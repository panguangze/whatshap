"""
Connectivity analysis:
want to implement a graph, where the nodes are the reads, here called element, which store the covered positions
and is identified, by a tuple of the minimal position ,maximal position and tht corresponding index.

The edges in the graph, include as weight the number of  positions covering both nodes.

In order to find the connected components in the graph, use ing either the BFS ir DFS to finde connected components,
but only with higher or even to the given weight for the edges.
During the search the components should be stored ina disjoint data structure.


"""
class Element:
    def __init__(self,positions,i ):
        self.positions = set(positions)
        self.min = min(positions)
        self.max=max(positions)
        self.index=i
        #include connections, or edges (Element,w) means edge to Element with weight w
        self.connections=[]
        #To which component this node is included
        self.component=None

    def add_connection(self,Element,w):
        #Connection=edge
        self.connections.append((Element,w))

    def get_positions(self):
        return self.positions

    def get_min(self):
        return self.min

    def get_max(self):
        return self.max

    def get_index(self):
        return self.index

    def get_connections(self):
        return self.connections

    def get_component(self):
        return  self.component

    def set_component(self,component):
        self.component=component

    #Difference between equal element and equal by sorting
    def equal_element(self,node):
        return  ( self.min==node.get_min() and
            self.max==node.get_max() and self.positions==node.get_positions()
                  and self.index==node.get_index())

    #Adding comparison methods for storting

    #=
    def __eq__(self, other):
        return( self.min==other.min and
            self.max==other.max and self.positions==other.positions)

    #x<y
    def __lt__(self, other):
        return (self.min<other.max or (self.min==other.min and self.max<self.max))

    #x<=y
    def __le__(self, other):
        return (self.min<self.other or (self.min==other.min and self.max<=self.max))

    #x!=y
    def __ne__(self, other):
        return (self.min==self.max and self.min==self.max and self.positions==other.positions)

    #x>y
    def __gt__(self, other):
        return (self.min>other.min or (self.min==other.min and self.max>other.max))

    #x>=y
    def __ge__(self, other):
        return (self.min>other.min or(self.min==other.min and self.max >=other.max))


class read_positions_graph:
    """
    Stores the graph by nodes and edges,
    also supports the DFS for finding connected components depending on the overlapping positions
    """

    def __init__(self, reads):
        self.nodes = []
        i=0
        while i != len(reads):
            positions = [variant.position for variant in reads[i]]
            new_node=Element(positions,i)
            self.nodes.append(new_node)
            i+=1


        #TODO do not know if neede , like in test case also without sorting the nodes are sorted
        sorted_nodes=sorted(self.nodes)
        sorted(self.nodes)


        #Stored the components of the graphs. for later use
        self.components=[]
        self.edges=[]
        #node which connections are added
        for node in self.nodes:
            node_min=node.get_min()
            node_max=node.get_max()
            for sec_node in self.nodes:
                #Do not need a self loop
                if sec_node.equal_element(node):
                    continue
                sec_min=sec_node.get_min()
                #TODO If sorted works could make here a cut if sec_node greater then node, or if node_max<sec_node_min ist
                #if sec_min>node_max:
                #	break
                sec_max=sec_node.get_max()
                #cases :
                    #1 : min & max are less than sec_min or bigger than sec_max => no intersection possible - no edges possible - printouts
                    #2 : min < sec_min  but max>=sec_min  - analysis overlap - independent of single or paired end
                    #3 : max>sec_max but min <sec_max - analysis overlap - independent of single or paired end
                    #4 :eiter s_min<min<s_max and and s_max>max>s_min or  min<s_min<max or max>s_max>min
                #Not implemented in this order
                #The nodes which include the same min and max, but different indices should have an edge
                if node_min<=sec_min:
                    #print('In first if')
                    if node_max>sec_max :
                        #sec_node positions are included into node positions
                        weight=len(sec_node.get_positions().intersection(node.get_positions()))
                        if ((sec_node,node,weight)not  in self.edges and weight!=0):
                            node.add_connection(sec_node,weight)
                            sec_node.add_connection(node,weight)
                            self.edges.append((node,sec_node,weight))
                    else:
                        if sec_min<=node_max:
                            #max of sec_node  is higher than node max, so  covering at least position node_max
                            weight=len(sec_node.get_positions().intersection(node.get_positions()))
                            if ((sec_node,node,weight)not  in self.edges and weight!=0):
                                sec_node.add_connection(node,weight)
                                node.add_connection(sec_node,weight)
                                self.edges.append((node,sec_node,weight))
                                #print('Added edge in the max of sec node highter than node max')
                                #print(node_min)
                                #print(node_max)
                                #print(sec_min)
                                #print(sec_max)
                if sec_min<node_min:
                    if node_max>sec_max:
                        weight=len(sec_node.get_positions().intersection(node.get_positions()))
                        if ((sec_node,node,weight)not  in self.edges and weight!=0):
                            node.add_connection(sec_node,weight)
                            sec_node.add_connection(node,weight)
                            self.edges.append((node,sec_node,weight))
                    else:
                        if sec_max>=node_min:
                            #node completely included into sec node
                            #Need to use at all stations the intersection because we could also llok at paired end reads, so the complete including of the read
                            #does not mean that alll positions of this read are also covered
                            weight=len(sec_node.get_positions().intersection(node.get_positions()))
                            if ((sec_node,node,weight)not  in self.edges and weight!=0):
                                node.add_connection(sec_node,weight)
                                sec_node.add_connection(node,weight)
                                self.edges.append((node,sec_node,weight))



    def get_edges(self):
        return self.edges

    def get_components(self):
        return self.components

    def get_nodes(self):
        return self.nodes

    def add_components(self,comp):
        #comp is an disjoint element
        self.components.append(comp)

class Connect_comp:
    #First each detected node is a connected comp, and then merged....
    #if there is no merge then the node remains in the not visited state and later analysed

    #Items = included_nodes: nodes in the connected components
    # blocks: dictionary of  keys (which are the degrees of the edges) and the values is the number of edges of this degree
    #           in the component
    #stored_for_later: Neighbour nodes, which in the first approach are not added, are stored there to later look again to them
    #Analyze_nodes: Nodes which need to be analyzed in the BFS approach
    #Length : Number of positions in the connected component
    #max = maximal position
    #min = min position




    def __init__(self,node,factor):
        #Initilaize the component, by the actual node and analyze depending on the factor if the neighbours need to be
        #stored for later or analyzed
        self.positions=node.get_positions() # set
        self.included_nodes =[node]
        self.stored_for_later=[]
        self.analyze_nodes=[]
        self.blocks={}
        neighbours=node.get_connections()

        for (element,w) in neighbours:
            if w>=factor:
                self.analyze_nodes.append(element)
                self.blocks[w]=1
            else:
                self.stored_for_later.append(element)
        #blocks store the number of connection between reads,depending on the degree w
        self.length= len(self.positions)
        self.max=max(self.positions)
        self.min=min(self.positions)


    def update_component(self,node):
        #node added to the actual component
        new_position=node.get_positions()
        #just union them
        self.positions=self.positions.union(new_position)
        self.included_nodes.append(node)
        self.length= len(self.positions)
        self.max=max(self.positions)
        self.min=min(self.positions)

    def expand_component(self,node,factor):
        #Expand the component by the neighbours of the given node
        neighbours=node.get_connections()
        for (element,w) in neighbours:
            if w>=factor:
                self.analyze_nodes.append(element)
                self.blocks[w]+=1
            else:
                self.stored_for_later.append(element)

    def get_length(self):
        return self.length

    def get_positions(self):
        return self.positions

    def get_included_nodes(self):
        return self.included_nodes

    def get_analyze_nodes(self):
        return self.analyze_nodes

    def get_blocks(self):
        return self.blocks

    def get_stored_for_later(self):
        return self.stored_for_later

    def get_min(self):
        return self.min

    def get_max(self):
        return self.max

#def BFS_search_in_read_graph(read_graph,connect_factor):
#    nodes=read_graph.get_nodes()
#    edges=read_graph.get_edges()


#def recursion_step(node,factor):
#    neighbours=node.get_connections()
#    component=node.get_component()
#    need_to_add=[(sec_node,weight) for (sec_node,weight) in neighbours if weight>=factor]
#    need_to_be_stored = [sec_node for (sec_node,weight) in neighbours if weight<factor]
#    component.add_stored(need_to_be_stored)
#    while len(need_to_add)!=0:
#        (connected_node,w)=need_to_add.pop()
#        connected_node.add_component(component)
#        component.update(connected_node.w)

def make_component(node,graph,factor):
    new_component=Connect_comp(node,factor)
    return new_component