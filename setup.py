"""
Following Cython's recommendations: If pre-generated .c extension files are
found, then Cython is not run, even if it is installed.
Pass --cython on the command line to force Cython to be run.

If .c files are not found, such as in a fresh Git checkout, then Cython is
always run.
"""
import sys
import os.path
from setuptools import setup, Extension
from distutils.version import LooseVersion
from whatshap import __version__

MIN_CYTHON_VERSION = '0.17'

if sys.version_info < (3, 3):
	sys.stdout.write("At least Python 3.3 is required.\n")
	sys.exit(1)


def out_of_date(extensions):
	"""
	Check whether any pyx source is newer than the corresponding generated
	C(++) source or whether any C(++) source is missing.
	"""
	for extension in extensions:
		for pyx in extension.sources:
			path, ext = os.path.splitext(pyx)
			if ext not in ('.pyx', '.py'):
				continue
			csource = path + ('.cpp' if extension.language == 'c++' else '.c')
			if not os.path.exists(csource) or (
				os.path.getmtime(pyx) > os.path.getmtime(csource)):
				return True
	return False


def no_cythonize(extensions, **_ignore):
	"""
	Change file extensions from .pyx to .c or .cpp.

	Copied from Cython documentation
	"""
	for extension in extensions:
		sources = []
		for sfile in extension.sources:
			path, ext = os.path.splitext(sfile)
			if ext in ('.pyx', '.py'):
				if extension.language == 'c++':
					ext = '.cpp'
				else:
					ext = '.c'
				sfile = path + ext
			sources.append(sfile)
		extension.sources[:] = sources
	return extensions


def cythonize_if_necessary(extensions):
	if '--cython' in sys.argv:
		sys.argv.remove('--cython')
	elif out_of_date(extensions):
		sys.stdout.write('At least one C source file is missing or out of date.\n')
	else:
		return no_cythonize(extensions)

	try:
		from Cython import __version__ as cyversion
	except ImportError:
		sys.stdout.write(
			"ERROR: Cython is not installed. Install at least Cython version " +
			str(MIN_CYTHON_VERSION) + " to continue.\n")
		sys.exit(1)
	if LooseVersion(cyversion) < LooseVersion(MIN_CYTHON_VERSION):
		sys.stdout.write(
			"ERROR: Your Cython is at version '" + str(cyversion) +
			"', but at least version " + str(MIN_CYTHON_VERSION) + " is required.\n")
		sys.exit(1)

	from Cython.Build import cythonize
	return cythonize(extensions)

extensions = [
	Extension('whatshap._core', sources=['whatshap/_core.pyx'], language='c++', extra_compile_args=["-std=c++0x"],),
]
extensions = cythonize_if_necessary(extensions)

setup(
	name = 'whatshap',
	version = __version__,
	author = '',
	author_email = '',
	url = 'https://bitbucket.org/whatshap/whatshap/',
	description = 'phase genomic variants using DNA sequencing reads',
	license = 'MIT',
	ext_modules = extensions,
	packages = ['whatshap'],
	#scripts = ['bin/...'],
	install_requires = ['pysam', 'PyVCF'],
	classifiers = [
		"Development Status :: 2 - Pre-Alpha",
		"Environment :: Console",
		"Intended Audience :: Science/Research",
		"License :: OSI Approved :: MIT License",
		"Natural Language :: English",
		"Programming Language :: Cython",
		"Programming Language :: Python",
		"Programming Language :: Python :: 3",
		"Topic :: Scientific/Engineering :: Bio-Informatics"
	]
)
