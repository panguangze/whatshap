import argparse


#!/usr/bin/env python3
"""
Phase variants in a VCF with the WhatsHap algorithm

Read a VCF and one or more files with phase information (BAM/CRAM or VCF phased
blocks) and phase the variants. The phased VCF is written to standard output.
"""
import logging
import sys
import platform

from argparse import SUPPRESS
from collections import defaultdict
from copy import deepcopy

from contextlib import ExitStack
from typing import Optional, List, TextIO, Union, Dict

from whatshap.vcf import VcfReader, PhasedVcfWriter, VcfError, VariantTable
from whatshap import __version__
from whatshap.core import (
    ReadSet,
    readselection,
    Pedigree,
    PedigreeDPTable,
    NumericSampleIds,
    PhredGenotypeLikelihoods,
)
from whatshap.graph import ComponentFinder
from whatshap.pedigree import (
    PedReader,
    mendelian_conflict,
    UniformRecombinationCostComputer,
    GeneticMapRecombinationCostComputer,
    find_recombination,
    ParseError,
    RecombinationCostComputer,
)
from whatshap.timer import StageTimer
from whatshap.utils import plural_s, warn_once
from whatshap.cli import CommandLineError, log_memory_usage, PhasedInputReader
from whatshap.merge import ReadMerger, DoNothingReadMerger, ReadMergerBase

__author__ = "Murray Patterson, Alexander Schönhuth, Tobias Marschall, Marcel Martin"

logger = logging.getLogger(__name__)


def find_components(phased_positions, reads, master_block=None, heterozygous_positions=None):
    """
    Return a dict that maps each variant position to the component it is in.
    Variants are considered to be in the same component if a read exists that
    covers both. A component is identified by the position of its leftmost
    variant.
    master_block -- List of positions in a "master block", i.e. all blocks containing
                    any of these positions are merged into one block.
    heterozygous_positions -- A dictionary mapping numeric sample ids to sets of
                              positions. Component building is then restricted to variants
                              at these positions. If none, all variants are used.
    """
    logger.debug("Finding connected components ...")
    assert phased_positions == sorted(phased_positions)

    # Find connected components.
    # A component is identified by the position of its leftmost variant.
    component_finder = ComponentFinder(phased_positions)
    phased_positions = set(phased_positions)
    for read in reads:
        if heterozygous_positions is None:
            positions = [
                variant.position for variant in read if variant.position in phased_positions
            ]
        else:
            positions = [
                variant.position
                for variant in read
                if (variant.position in phased_positions)
                and (variant.position in heterozygous_positions[read.sample_id])
            ]
        for position in positions[1:]:
            component_finder.merge(positions[0], position)
    if master_block is not None:
        for position in master_block[1:]:
            component_finder.merge(master_block[0], position)
    components = {position: component_finder.find(position) for position in phased_positions}
    return components


def find_largest_component(components):
    """
    Determine the largest component and return a sorted list of positions
    contained in it.
    components -- dictionary mapping positin to block_id as returned by find_components.
    """
    blocks = defaultdict(list)
    for position, block_id in components.items():
        blocks[block_id].append(position)
    largest = []
    for block in blocks.values():
        if len(block) > len(largest):
            largest = block
    largest.sort()
    return largest


def best_case_blocks(reads):
    """
    Given a list of core reads, determine the number of phased blocks that
    would result if each variant were actually phased.

    Return the number of connected components and non-singleton components.
    """
    positions = set()
    for read in reads:
        for variant in read:
            positions.add(variant.position)
    component_finder = ComponentFinder(positions)
    for read in reads:
        read_positions = [variant.position for variant in read]
        for position in read_positions[1:]:
            component_finder.merge(read_positions[0], position)
    # A dict that maps each component to the number of variants it contains
    component_sizes = defaultdict(int)
    for position in positions:
        component_sizes[component_finder.find(position)] += 1
    non_singletons = [component for component, size in component_sizes.items() if size > 1]
    return len(component_sizes), len(non_singletons)


def select_reads(readset, max_coverage, preferred_source_ids):
    logger.info(
        "Reducing coverage to at most %dX by selecting most informative reads ...", max_coverage
    )
    selected_indices = readselection(readset, max_coverage, preferred_source_ids)
    selected_reads = readset.subset(selected_indices)
    logger.info(
        "Selected %d reads covering %d variants",
        len(selected_reads),
        len(selected_reads.get_positions()),
    )

    return selected_reads

def setup_pedigree(ped_path, samples):
    """
    Read in PED file to set up list of relationships.

    Return a pair (trios, pedigree_samples), where trios is a list of Trio
    objects and pedigree_samples is the set of all samples that are mentioned
    in the PED file (as individual, mother or father).

    ped_path -- path to PED file
    samples -- samples to be phased
    """
    trios = []
    pedigree_samples = set()
    for trio in PedReader(ped_path):
        if trio.child is None or trio.mother is None or trio.father is None:
            warn_once(
                logger,
                "Relationship %s/%s/%s ignored because at least one of the individuals is unknown.",
                trio.child,
                trio.mother,
                trio.father,
            )
            continue
        # if at least one individual is not in samples, skip trio
        if (
            (trio.mother not in samples)
            or (trio.father not in samples)
            or (trio.child not in samples)
        ):
            # happens in case --ped and --samples are used
            warn_once(
                logger,
                "Relationship %s/%s/%s ignored because at least one of the "
                "individuals was not given by --samples.",
                trio.child,
                trio.mother,
                trio.father,
            )
            continue

        trios.append(trio)
        pedigree_samples.add(trio.child)
        pedigree_samples.add(trio.father)
        pedigree_samples.add(trio.mother)

    return trios, pedigree_samples

def phase(sample: str, variant_table: VariantTable, chromosome, families, family_trios, max_coverage, timers: StageTimer):
    chromosome = variant_table.chromosome
    # These two variables hold the phasing results for all samples
    superreads, components = dict(), dict()

    # Iterate over all families to process, i.e. a separate DP table is created
    # for each family.
    # TODO: Can the body of this loop be factored out into a phase_family function?
    for representative_sample, family in sorted(families.items()):
        logger.info("---- Processing family with individuals: %s", ",".join(family))
        max_coverage_per_sample = max(1, max_coverage // len(family))
        logger.info("Using maximum coverage per sample of %dX", max_coverage_per_sample)
        trios = family_trios[representative_sample]
        assert len(family) == 1 or len(trios) > 0

        homozygous_positions, phasable_variant_table = find_phaseable_variants(
            family, trios, variant_table
        )

        # Get the reads belonging to each sample
        readsets = dict()  # TODO this could become a list
        for sample in family:
            with timers("read_bam"):
                readset, vcf_source_ids = phased_input_reader.read(
                    chromosome, phasable_variant_table.variants, sample
                )

            # TODO: Read selection done w.r.t. all variants, where using heterozygous
            #  variants only would probably give better results.
            with timers("select"):
                readset = readset.subset(
                    [i for i, read in enumerate(readset) if len(read) >= 2]
                )
                logger.info(
                    "Kept %d reads that cover at least two variants each", len(readset)
                )
                merged_reads = read_merger.merge(readset)
                selected_reads = select_reads(
                    merged_reads,
                    max_coverage_per_sample,
                    preferred_source_ids=vcf_source_ids,
                )

            readsets[sample] = selected_reads

        all_reads = merge_readsets(readsets)

        # Determine which variants can (in principle) be phased
        accessible_positions = sorted(all_reads.get_positions())
        logger.info(
            "Variants covered by at least one phase-informative "
            "read in at least one individual after read selection: %d",
            len(accessible_positions),
        )
        if len(family) > 1 and genetic_haplotyping:
            # In case of genetic haplotyping, also retain all positions homozygous
            # in at least one individual (because they might be phased based on genotypes)
            accessible_positions = sorted(
                set(accessible_positions).union(homozygous_positions)
            )
            logger.info(
                "Variants either covered by phase-informative read or homozygous "
                "in at least one individual: %d",
                len(accessible_positions),
            )

        # Keep only accessible positions
        phasable_variant_table.subset_rows_by_position(accessible_positions)
        assert len(phasable_variant_table.variants) == len(accessible_positions)

        pedigree = create_pedigree(
            default_gq,
            family,
            numeric_sample_ids,
            phasable_variant_table,
            trios,
        )
        # recombination_costs = recombination_cost_computer.compute(accessible_positions)

        # Finally, run phasing algorithm
        # with timers("phase"):
        #     problem_name = "MEC" if len(family) == 1 else "PedMEC"
        #     logger.info(
        #         "Phasing %d sample%s by solving the %s problem ...",
        #         len(family),
        #         plural_s(len(family)),
        #         problem_name,
        #     )

        #     dp_table: Union[HapChatCore, PedigreeDPTable]

            # dp_table = PedigreeDPTable(
            #     all_reads,
            #     recombination_costs,
            #     pedigree,
            #     accessible_positions,
            # )

        #     superreads_list, transmission_vector = dp_table.get_super_reads()
        #     logger.info("%s cost: %d", problem_name, dp_table.get_optimal_cost())

        # with timers("components"):
        #     overall_components = compute_overall_components(
        #         accessible_positions,
        #         all_reads,
        #         family,
        #         genetic_haplotyping,
        #         homozygous_positions,
        #         numeric_sample_ids,
        #         superreads_list,
        #     )
        #     log_component_stats(overall_components, len(accessible_positions))
        # Superreads in superreads_list are in the same order as individuals were added to the pedigree
        # for sample, sample_superreads in zip(family, superreads_list):
        #     superreads[sample] = sample_superreads
        #     assert len(sample_superreads) == 2
        #     assert (
        #         sample_superreads[0].sample_id
        #         == sample_superreads[1].sample_id
        #         == numeric_sample_ids[sample]
        #     )
        #     # identical for all samples
        #     components[sample] = overall_components

    logger.debug("Chromosome %r finished", chromosome)


def run_whatshap(
    variant_file: str,
    output: TextIO = sys.stdout,
    chromosomes: Optional[List[str]] = None,
    indels: bool = True,
    mapping_quality: int = 20,
    max_coverage: int = 15,
    ped: Optional[str] = None,
    recombrate: float = 1.26,
    genmap: Optional[str] = None,
    genetic_haplotyping: bool = True,
    default_gq: int = 30,
    use_ped_samples: bool = False,
    tag: str = "PS",
):
    """
    Run WhatsHap.

    phase_input_files -- list of paths to BAM/CRAM/VCF files
    variant_file -- path to input VCF
    reference -- path to reference FASTA. If False: skip realignment. If None: complain if reference needed.
    output -- path to output VCF or a file-like object
    samples -- names of samples to phase. an empty list means: phase all samples
    chromosomes -- names of chromosomes to phase. an empty list means: phase all chromosomes
    ignore_read_groups
    mapping_quality -- discard reads below this mapping quality
    read_merging -- whether or not to merge reads
    read_merging_error_rate -- probability that a nucleotide is wrong
    read_merging_max_error_rate -- max error rate on edge of merge graph considered
    read_merging_positive_threshold -- threshold on the ratio of the two probabilities
    read_merging_negative_threshold -- threshold on the opposite ratio of positive threshold
    max_coverage
    distrust_genotypes
    include_homozygous
    genetic_haplotyping -- in ped mode, merge disconnected blocks based on genotype status
    recombination_list_filename -- filename to write putative recombination events to
    tag -- How to store phasing info in the VCF, can be 'PS' or 'HP'
    algorithm -- algorithm to use, can be 'whatshap' or 'hapchat'
    gl_regularizer -- float to be passed as regularization constant to GenotypeLikelihoods.as_phred
    gtchange_list_filename -- filename to write list of changed genotypes to
    default_gq -- genotype likelihood to be used when GL or PL not available
    write_command_line_header -- whether to add a ##commandline header to the output VCF
    """
    timers = StageTimer()
    logger.info(f"This is WhatsHap {__version__} running under Python {platform.python_version()}")
    numeric_sample_ids = NumericSampleIds()
    command_line: Optional[str]
    command_line = None
    read_merger = DoNothingReadMerger()
    with ExitStack() as stack:
        try:
            vcf_writer = stack.enter_context(
                PhasedVcfWriter(
                    command_line=command_line,
                    in_path=variant_file,
                    out_file=output,
                    tag=tag,
                    indels=indels,
                )
            )
        except (OSError, VcfError) as e:
            raise CommandLineError(e)

        phased_input_reader = stack.enter_context(
            PhasedInputReader(
                variant_file,
                numeric_sample_ids,
                mapq_threshold=mapping_quality,
                indels=indels,
            )
        )
        show_phase_vcfs = phased_input_reader.has_vcfs

        # Only read genotype likelihoods from VCFs when distrusting genotypes
        # vcf_reader = stack.enter_context(
        #     VcfReader(variant_file, indels=indels)
        # )
        # samples = vcf_reader.samples

        # if --use-ped-samples is set, use only samples from PED file
        if ped and use_ped_samples:
            samples = PedReader(ped).samples()

        # raise_if_any_sample_not_in_vcf(vcf_reader, samples)

        # recombination_cost_computer = make_recombination_cost_computer(ped, genmap, recombrate)

        families, family_trios = setup_families(samples, ped, max_coverage)
        del samples
        for trios in family_trios.values():
            for trio in trios:
                # Ensure that all mentioned individuals have a numeric id
                _ = numeric_sample_ids[trio.child]

        read_list = None

        with timers("parse_phasing_vcfs"):
            # TODO should this be done in PhasedInputReader.__init__?
            phased_input_reader.read_vcfs()

        superreads: Dict[str, ReadSet]
        components: Dict
        for variant_table in timers.iterate("parse_vcf", vcf_reader):
            chromosome = variant_table.chromosome
            if (not chromosomes) or (chromosome in chromosomes):
                logger.info("======== Working on chromosome %r", chromosome)
            else:
                logger.info(
                    "Leaving chromosome %r unchanged (present in VCF but not requested by option --chromosome)",
                    chromosome,
                )
                with timers("write_vcf"):
                    superreads, components = dict(), dict()
                    vcf_writer.write(chromosome, superreads, components)
                continue

            # These two variables hold the phasing results for all samples
            superreads, components = dict(), dict()

            # Iterate over all families to process, i.e. a separate DP table is created
            # for each family.
            # TODO: Can the body of this loop be factored out into a phase_family function?
            for representative_sample, family in sorted(families.items()):
                if len(family) == 1:
                    logger.info("---- Processing individual %s", representative_sample)
                else:
                    logger.info("---- Processing family with individuals: %s", ",".join(family))
                max_coverage_per_sample = max(1, max_coverage // len(family))
                logger.info("Using maximum coverage per sample of %dX", max_coverage_per_sample)
                trios = family_trios[representative_sample]
                assert len(family) == 1 or len(trios) > 0

                homozygous_positions, phasable_variant_table = find_phaseable_variants(
                    family, trios, variant_table
                )

                # Get the reads belonging to each sample
                readsets = dict()  # TODO this could become a list
                for sample in family:
                    with timers("read_bam"):
                        readset, vcf_source_ids = phased_input_reader.read(
                            chromosome, phasable_variant_table.variants, sample
                        )

                    # TODO: Read selection done w.r.t. all variants, where using heterozygous
                    #  variants only would probably give better results.
                    with timers("select"):
                        readset = readset.subset(
                            [i for i, read in enumerate(readset) if len(read) >= 2]
                        )
                        logger.info(
                            "Kept %d reads that cover at least two variants each", len(readset)
                        )
                        merged_reads = read_merger.merge(readset)
                        selected_reads = select_reads(
                            merged_reads,
                            max_coverage_per_sample,
                            preferred_source_ids=vcf_source_ids,
                        )

                    readsets[sample] = selected_reads

                all_reads = merge_readsets(readsets)

                # Determine which variants can (in principle) be phased
                accessible_positions = sorted(all_reads.get_positions())
                logger.info(
                    "Variants covered by at least one phase-informative "
                    "read in at least one individual after read selection: %d",
                    len(accessible_positions),
                )
                if len(family) > 1 and genetic_haplotyping:
                    # In case of genetic haplotyping, also retain all positions homozygous
                    # in at least one individual (because they might be phased based on genotypes)
                    accessible_positions = sorted(
                        set(accessible_positions).union(homozygous_positions)
                    )
                    logger.info(
                        "Variants either covered by phase-informative read or homozygous "
                        "in at least one individual: %d",
                        len(accessible_positions),
                    )

                # Keep only accessible positions
                phasable_variant_table.subset_rows_by_position(accessible_positions)
                assert len(phasable_variant_table.variants) == len(accessible_positions)

                pedigree = create_pedigree(
                    default_gq,
                    family,
                    numeric_sample_ids,
                    phasable_variant_table,
                    trios,
                )
                # recombination_costs = recombination_cost_computer.compute(accessible_positions)

                # Finally, run phasing algorithm
                # with timers("phase"):
                #     problem_name = "MEC" if len(family) == 1 else "PedMEC"
                #     logger.info(
                #         "Phasing %d sample%s by solving the %s problem ...",
                #         len(family),
                #         plural_s(len(family)),
                #         problem_name,
                #     )

                #     dp_table: Union[HapChatCore, PedigreeDPTable]

                    # dp_table = PedigreeDPTable(
                    #     all_reads,
                    #     recombination_costs,
                    #     pedigree,
                    #     accessible_positions,
                    # )

                #     superreads_list, transmission_vector = dp_table.get_super_reads()
                #     logger.info("%s cost: %d", problem_name, dp_table.get_optimal_cost())

                # with timers("components"):
                #     overall_components = compute_overall_components(
                #         accessible_positions,
                #         all_reads,
                #         family,
                #         genetic_haplotyping,
                #         homozygous_positions,
                #         numeric_sample_ids,
                #         superreads_list,
                #     )
                #     log_component_stats(overall_components, len(accessible_positions))
                # Superreads in superreads_list are in the same order as individuals were added to the pedigree
                # for sample, sample_superreads in zip(family, superreads_list):
                #     superreads[sample] = sample_superreads
                #     assert len(sample_superreads) == 2
                #     assert (
                #         sample_superreads[0].sample_id
                #         == sample_superreads[1].sample_id
                #         == numeric_sample_ids[sample]
                #     )
                #     # identical for all samples
                #     components[sample] = overall_components

            logger.debug("Chromosome %r finished", chromosome)

    log_time_and_memory_usage(timers, show_phase_vcfs=show_phase_vcfs)


def compute_overall_components(
    accessible_positions,
    all_reads,
    family,
    genetic_haplotyping,
    homozygous_positions,
    numeric_sample_ids,
    superreads_list,
):
    master_block = None
    heterozygous_positions_by_sample = None
    # If we distrusted genotypes, we need to re-determine which sites are homo-/heterozygous after phasing
    if len(family) > 1 and genetic_haplotyping:
        master_block = sorted(set(homozygous_positions).intersection(set(accessible_positions)))
    return find_components(
        accessible_positions, all_reads, master_block, heterozygous_positions_by_sample
    )


def log_component_stats(components, n_accessible_positions):
    n_phased_blocks = len(set(components.values()))
    logger.info(f"No. of phased blocks: {n_phased_blocks}")
    largest = find_largest_component(components)
    if not largest:
        return
    logger.info(
        f"Largest block contains {len(largest)} variants"
        f" ({len(largest) / n_accessible_positions:.1%} of accessible variants)"
        f" between position {largest[0] + 1} and {largest[-1] + 1}"
    )


def log_best_case_phasing_info(readset, selected_reads):
    (n_best_case_blocks, n_best_case_nonsingleton_blocks) = best_case_blocks(readset)
    (n_best_case_blocks_cov, n_best_case_nonsingleton_blocks_cov) = best_case_blocks(selected_reads)
    logger.info(
        "Best-case phasing would result in %d non-singleton phased blocks (%d in total)",
        n_best_case_nonsingleton_blocks,
        n_best_case_blocks,
    )
    logger.info(
        "... after read selection: %d non-singleton phased blocks (%d in total)",
        n_best_case_nonsingleton_blocks_cov,
        n_best_case_blocks_cov,
    )


def raise_if_any_sample_not_in_vcf(vcf_reader, samples):
    vcf_sample_set = set(vcf_reader.samples)
    for sample in samples:
        if sample not in vcf_sample_set:
            raise CommandLineError(
                "Sample {!r} requested on command-line not found in VCF".format(sample)
            )


def setup_families(samples, ped, max_coverage):
    """
    Return families, family_trios pair.

    families maps a family representative to a list of family members

    family_trios maps a family representative to a list of trios in this family
    """

    # list of all trios across all families
    all_trios = dict()

    # Keep track of connected components (aka families) in the pedigree
    family_finder = ComponentFinder(samples)

    if ped:
        all_trios, pedigree_samples = setup_pedigree(ped, samples)
        for trio in all_trios:
            family_finder.merge(trio.father, trio.child)
            family_finder.merge(trio.mother, trio.child)

    # map family representatives to lists of family members
    families = defaultdict(list)
    for sample in samples:
        families[family_finder.find(sample)].append(sample)

    # map family representatives to lists of trios for this family
    family_trios = defaultdict(list)
    for trio in all_trios:
        family_trios[family_finder.find(trio.child)].append(trio)
    logger.info(
        "Working on %d%s samples from %d famil%s",
        len(samples),
        plural_s(len(samples)),
        len(families),
        "y" if len(families) == 1 else "ies",
    )

    largest_trio_count = max([0] + [len(trio_list) for trio_list in family_trios.values()])
    if max_coverage + 2 * largest_trio_count > 23:
        logger.warning(
            "The maximum coverage is too high! "
            "WhatsHap may take a long time to finish and require a huge amount of memory."
        )
    return families, family_trios


def make_recombination_cost_computer(
    ped: Optional[str], genmap: Optional[str], recombrate: float
) -> RecombinationCostComputer:
    if ped and genmap:
        logger.info("Using region-specific recombination rates from genetic map %s.", genmap)
        try:
            return GeneticMapRecombinationCostComputer(genmap)
        except ParseError as e:
            raise CommandLineError(e)
    else:
        if ped:
            logger.info("Using uniform recombination rate of %g cM/Mb.", recombrate)
        return UniformRecombinationCostComputer(recombrate)


def find_phaseable_variants(family, trios, variant_table: VariantTable):
    # variant indices with at least one missing genotype
    missing_genotypes = set()
    # variant indices with at least one heterozygous genotype
    heterozygous = set()
    # variant indices with at least one homozygous genotype
    homozygous = set()
    # determine which variants have missing/heterozygous/homozygous genotypes in any sample
    for sample in family:
        genotypes = variant_table.genotypes_of(sample)
        for index, gt in enumerate(genotypes):
            if gt.is_none():
                missing_genotypes.add(index)
            elif not gt.is_homozygous():
                heterozygous.add(index)
            else:
                assert gt.is_diploid_and_biallelic()
                homozygous.add(index)
    # determine which variants have Mendelian conflicts
    # variant indices with at least one Mendelian conflict
    mendelian_conflicts = find_mendelian_conflicts(trios, variant_table)
    # retain variants that are heterozygous in at least one individual (anywhere in the pedigree)
    # and do not have neither missing genotypes nor Mendelian conflict
    to_retain = heterozygous
    to_retain = to_retain.difference(missing_genotypes).difference(mendelian_conflicts)
    # discard every variant that is not to be retained
    to_discard = set(range(len(variant_table))).difference(to_retain)
    # Determine positions of selected variants that are homozygous in at least one individual.
    # These are used later to merge blocks containing these variants into one block (since
    # the are connected by "genetic haplotyping").
    homozygous_positions = [
        variant_table.variants[i].position for i in to_retain.intersection(homozygous)
    ]
    phasable_variant_table = deepcopy(variant_table)
    # Remove calls to be discarded from variant table
    phasable_variant_table.remove_rows_by_index(to_discard)
    logger.info("Number of variants skipped due to missing genotypes: %d", len(missing_genotypes))
    if len(family) == 1:
        logger.info(
            "Number of remaining%s variants: %d",
            "",
            len(phasable_variant_table),
        )
    else:
        logger.info(
            "Number of variants skipped due to Mendelian conflicts: %d", len(mendelian_conflicts)
        )
        logger.info(
            "Number of remaining variants heterozygous in at least one individual: %d",
            len(phasable_variant_table),
        )
    return homozygous_positions, phasable_variant_table


def log_time_and_memory_usage(timers, show_phase_vcfs):
    total_time = timers.total()
    logger.info("\n== SUMMARY ==")
    log_memory_usage()
    # fmt: off
    logger.info("Time spent reading BAM/CRAM:                 %6.1f s", timers.elapsed("read_bam"))
    logger.info("Time spent parsing VCF:                      %6.1f s", timers.elapsed("parse_vcf"))
    if show_phase_vcfs:
        logger.info("Time spent parsing input phasings from VCFs: %6.1f s", timers.elapsed("parse_phasing_vcfs"))
    logger.info("Time spent selecting reads:                  %6.1f s", timers.elapsed("select"))
    logger.info("Time spent phasing:                          %6.1f s", timers.elapsed("phase"))
    logger.info("Time spent writing VCF:                      %6.1f s", timers.elapsed("write_vcf"))
    logger.info("Time spent finding components:               %6.1f s", timers.elapsed("components"))
    logger.info("Time spent on rest:                          %6.1f s", total_time - timers.sum())
    logger.info("Total elapsed time:                          %6.1f s", total_time)
    # fmt: on


def merge_readsets(readsets) -> ReadSet:
    all_reads = ReadSet()
    for sample, readset in readsets.items():
        for read in readset:
            assert read.is_sorted(), "Add a read.sort() here"
            all_reads.add(read)
    all_reads.sort()
    return all_reads


def create_pedigree(
    default_gq,
    family,
    numeric_sample_ids,
    phasable_variant_table,
    trios,
):
    pedigree = Pedigree(numeric_sample_ids)
    for sample in family:
        # If distrusting genotypes, we pass genotype likelihoods on to pedigree object

        genotype_likelihoods = None
        pedigree.add_individual(
            sample, phasable_variant_table.genotypes_of(sample), genotype_likelihoods
        )
    for trio in trios:
        pedigree.add_relationship(father_id=trio.father, mother_id=trio.mother, child_id=trio.child)
    return pedigree


def find_mendelian_conflicts(trios, variant_table):
    mendelian_conflicts = set()
    for trio in trios:
        genotypes_mother = variant_table.genotypes_of(trio.mother)
        genotypes_father = variant_table.genotypes_of(trio.father)
        genotypes_child = variant_table.genotypes_of(trio.child)

        for index, (gt_mother, gt_father, gt_child) in enumerate(
            zip(genotypes_mother, genotypes_father, genotypes_child)
        ):
            if (not gt_mother.is_none()) and (not gt_father.is_none()) and (not gt_child.is_none()):
                if mendelian_conflict(gt_mother, gt_father, gt_child):
                    mendelian_conflicts.add(index)
    return mendelian_conflicts

def main():
    # phase_input_files: List[str],
    # variant_file: str,
    # output: TextIO = sys.stdout,
    # chromosomes: Optional[List[str]] = None,
    # indels: bool = True,
    # mapping_quality: int = 20,
    # max_coverage: int = 15,
    # ped: Optional[str] = None,
    # recombrate: float = 1.26,
    # genmap: Optional[str] = None,
    # genetic_haplotyping: bool = True,
    # default_gq: int = 30,
    # use_ped_samples: bool = False,
    # tag: str = "PS",
    parser = argparse.ArgumentParser("trio phase")
    parser.add_argument(
        '-v', help='merged VCF file', required=True, dest='vcf_file')
    parser.add_argument(
        '-p', help='pedigree file', required=True, dest='ped_fle')
    parser.add_argument(
        '-o', help='out phased vcf file', required=True, dest='out_file')
    args = parser.parse_args()

    run_whatshap(args.vcf_file, args.ped_fle, args. out_file)


if __name__ == "__main__":
    main()

