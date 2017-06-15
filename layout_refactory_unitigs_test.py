import networkx as nx
import sys
import operator
import pickle
import argparse
import os

parser = argparse.ArgumentParser()
#parser.add_argument('-a','--assembly', help='Contig assembly', required=False)
parser.add_argument('-x','--graph', help='GFA file for assembly graph', required=False,default='abc')
parser.add_argument('-l','--links', help='Links sorted by relative score', required=True)
parser.add_argument('-c','--cutoff', help='Minimum length contig to consider for scaffolding', required=False)
parser.add_argument('-i','--iteration',help='Iteration number',required=False,default=1)
parser.add_argument('-u','--unitigs',help='Bed file for unitig to contig tiline',required=False,default = 'abc')
parser.add_argument('-t','--tenx',help='Links obtained from 10x file',required=False,default = 'abc')
parser.add_argument('-d','--directory',help='Output Directory',required=False,default='out')
args = parser.parse_args()


G = nx.Graph()
contigs = set()
all_G = nx.Graph()
contig_length = {}

contig_lengths_original = {}
with open(args.directory+'/scaffold_length_iteration_1','r') as f:
    for line in f:
        attrs = line.split()
        contig_length[attrs[0]] = int(attrs[1])

OVl_G = nx.Graph()

'''
Given two nodes, it finds the ratio of shortest with second shortest path
'''
def get_best_path(start, end):
   ratio=0
   minpath=-1
   ori=""

   for link in ["B:E", "B:B", "E:B", "E:E"]:
      edgeType=link.split(":")
      p = get_path_two(start+"-"+edgeType[0], end+"-"+edgeType[1])
      if minpath == -1 or len(p) < minpath:
         if (minpath != -1):
            ratio = (float)(len(p))/minpath
         minpath=len(p)
         ori=link
   return [minpath, ori, ratio]

'''
Helper function for get_best_path(start,end)
'''
def get_path_two(start,end):
    v1 = start.split('-')
    v2 = end.split('-')

    if v1[1] == 'B':
        or1 = 'REV'
    else:
        or1 = 'FOW'

    if v2[1] == 'B':
        or2 = 'FOW'
    else:
        or2 = 'REV'
    start = v1[0]+'$'+or1
    end = v2[0] + '$'+or2
    paths = nx.single_source_dijkstra(OVl_G,source=start,target=end)
  
    paths = paths[1]
    path = []
    if end in paths:
        path = paths[end]
    else:
        return "NO PATH FOUND"
    ret = []
    for node in path:
        n = node.split('$')
        if n[1] == 'REV':
            ret.append(n[0]+'-E')
            ret.append(n[0]+'-B')
        else:
            ret.append(n[0]+'-B')
            ret.append(n[0]+'-E')
    return ret

'''
Loads GFA file for either unitigs or contigs, given the file path
'''
def load_GFA(path):
    #f = open(path,'r')
    with open(path,'r') as f:
        for line in f:
           line = line.strip().split()

           #S   tig00000001 *   LN:i:20388
           #L   tig00023281 -   tig00008904 +   21556M  cv:A:F
           if (line[0] == "S"):
              OVl_G.add_node(line[1]+"$REV", name=line[1], length=line[3].split(":")[2])
              OVl_G.add_node(line[1]+"$FOW", name=line[1], length=line[3].split(":")[2])
           if (line[0] == "L"):
              #BUG BUG BUG in GFA inverting edges on input
              # ctg1=line[3]+"_pilon"
              # ctg2=line[1]+"_pilon"
              # tmp=line[2]
              # line[2]=line[4]
              # line[4]=tmp
              #assume orientations are correct
              ctg1 = line[1]
              ctg2 = line[3]
              if (line[2] == "+" and line[4] == "+"):
                OVl_G.add_edge(ctg1+"$FOW", ctg2+"$FOW", weight=1)
                OVl_G.add_edge(ctg2+"$REV", ctg1+"$REV", weight=1)
              if (line[2] == "+" and line[4] == "-"):
                OVl_G.add_edge(ctg1+"$FOW", ctg2+"$REV", weight=1)
                OVl_G.add_edge(ctg2+"$FOW", ctg1+"$REV", weight=1)
              if (line[2] == "-" and line[4] == "+"):
                OVl_G.add_edge(ctg1+"$REV", ctg2+"$FOW", weight=1)
                OVl_G.add_edge(ctg2+"$REV", ctg1+"$FOW", weight=1)
              if (line[2] == "-" and line[4] == "-"):
                OVl_G.add_edge(ctg1+"$REV", ctg2+"$REV", weight=1)
                OVl_G.add_edge(ctg2+"$FOW", ctg1+"$FOW", weight=1)

'''
Find reverse complement of current orientation of contig
'''
def reverse_complement(contig):
    c1, o1 = contig.split(':')
    if o1 == 'B':
        return c1+':E'
    else:
        return c1+':B'


'''
Creates a graph based on unitig tiling bed file. 
Currently, all the orientations in the bed file are forward 
'''
def load_unitig_mapping():
    unitig_graph = nx.DiGraph()
    prev_line = ''
    with open(args.unitigs,'r') as f:
        for line in f:
            if prev_line == '':
                prev_line = line
                continue
            attrs = line.split()
            prev_attrs = prev_line.split()
            contig = attrs[0]
            prev_contig = prev_attrs[0]
            if prev_contig == contig:
                #this is to take into account the naming 
                prev_unitig = 'tig'+str(prev_attrs[3][3:])
                curr_unitig = 'tig'+str(attrs[3][3:])
                prev_orient = prev_attrs[-1]
                curr_orient = attrs[-1]
                if prev_orient == '+' and curr_orient == '+':
                    unitig_graph.add_edge(prev_unitig+':E',curr_unitig+':B')
                    unitig_graph.add_edge(curr_unitig+':E', prev_unitig+':B')
                if prev_orient == '+' and curr_orient == '-':
                    unitig_graph.add_edge(prev_unitig+':E',curr_unitig+':E')
                    unitig_graph.add_edge(curr_unitig+':B',prev_unitig+':B')
                if prev_orient == '-' and curr_unitig == '-':
                    unitig_graph.add_edge(prev_unitig+':B',curr_unitig+':E')
                    unitig_graph.add_edge(curr_unitig+':B', prev_unitig+':E')
                if prev_orient == '-' and curr_orient == '+':
                    unitig_graph.add_edge(prev_unitig+':B',curr_unitig+':B')
                    unitig_graph.add_edge(curr_unitig+':E',prev_unitig+':E')

            prev_line = line
    print >> sys.stderr, 'Unitig tiling graph loaded, nodes = ' + str(len(unitig_graph.nodes())) + ' edges = ' + str(len(unitig_graph.edges()))
    return unitig_graph


'''
Loads the 10x graph based on the bed links. Note that in this graph, only best weight edge between 
nodes is stored. Cutoff is for the score above which we want to keep the links
'''
def load_tenx_graph(cutoff):
    G_tenx = nx.Graph()
    with open(args.tenx,'r') as f:
        for line in f:
            attrs = line.split()
            v1 = attrs[0]
            v2 = attrs[1]
            if float(attrs[4]) >= cutoff and v1 not in tenx_graph.nodes() and v2 not in tenx_graph.nodes():
                G_tenx.add_edge(v1,v2,score=float(attrs[4]))

    return G_tenx


iteration = int(args.iteration)

'''
If GFA is given as an argument, load it
'''
shortest_paths = {}
if args.graph !='abc':
    print >> sys.stderr, 'started loading GFA'
    load_GFA(args.graph)
    #shortest_paths = nx.all_pairs_dijkstra_path_length(OVl_G)
print >> sys.stderr, 'Finished loading GFA, nodes =  ' + str(len(OVl_G.nodes())) + ' edges = ' + str(len(OVl_G.edges()))


'''
Load the contigs/scaffolds from the previous iteration if possible. 
Also store the mapping for the contigs to the end of the scaffolds to scaffold id
'''
previous_scaffolds = {}
contig2scaffold = {}
if int(args.iteration) > 1:
    try:
        previous_scaffolds = pickle.load(open(args.directory+'/scaffolds_iteration_'+str(iteration-1)+'.p','r'))
        for key in previous_scaffolds:
            contigs = previous_scaffolds[key]
            first = contigs[0].split(':')[0]
            last = contigs[-1].split(':')[0]
            contig2scaffold[first] = key
            contig2scaffold[last] = key
    except:
        previous_scaffolds = {}

'''
Load the unitig tiling and 10x graphs if possible.
'''
unitig_graph = nx.DiGraph()
tenx_graph = nx.Graph()

if args.unitigs != 'abc':
    unitig_graph = load_unitig_mapping()

'''
Load tenx graph this way only if it is a higher iteration 
'''
if args.tenx != 'abc' and iteration > 1:
    tenx_graph = load_tenx_graph(50)



#Load the best score graph first, keep log of used and unused links
print >> sys.stderr,  'Loading Hi-C links '
'''
Counts to keep log of each types of edge loaded in the graph
'''
hic_edges = 0
tiling_edges = 0
gfa_edges = 0
tenx_links = 0


'''
Given two scaffolds and either a unitig tiling or 10x graph, it finds out 
which orientations two scaffolds should be there if possible from the graph. 
Only applicable after first iteration.
'''

def test_edge(scaffold_first,scaffold_second,G_test,type):
    first_first = scaffold_first[0]
    last_first = scaffold_first[-1]
    first_second = scaffold_second[0]
    last_second = scaffold_second[-1]
    
    if G_test.has_edge(last_first,first_second):
        v1 = c1+':E'
        v2 = c2 +':B'
        if v1 not in G.nodes() and v2 not in G.nodes():
            G.add_edge(v1,v2,score=2,linktype=type)
            contigs.add(c1)
            contigs.add(c2)

    if G_test.has_edge(reverse_complement(first_first),first_second):
        v1 = c1 + ':B'
        v2 = c2 + ':B'
        if v1 not in G.nodes() and v2 not in G.nodes():
            G.add_edge(v1,v2,score=2,linktype=type)
            contigs.add(c1)
            contigs.add(c2)
            return 

    if G_test.has_edge(last_first,reverse_complement(last_second)):
        v1 = c1 + ':E'
        v2 = c2 + ':E'
        if v1 not in G.nodes() and v2 not in G.nodes():
            G.add_edge(v1,v2,score=2,linktype=type)
            contigs.add(c1)
            contigs.add(c2)
            return

    if G_test.has_edge(reverse_complement(first_first),reverse_complement(last_second)):
        v1 = c1 + ':E'
        v2 = c2 + ':B'
        if v1 not in G.nodes() and v2 not in G.nodes():
            G.add_edge(v1,v2,score=2,linktype=type)
            contigs.add(c1)
            contigs.add(c2)
            return

            #tiling_edges += 1


'''
This function loads 10x links 
'''

def load_tenx_links():
    if iteration == 1:
        if args.tenx != 'abc':
            with open(args.tenx,'r') as f:
                for line in f:
                    attrs = line.split()
                    v1 = attrs[0] 
                    v2 = attrs[1]
                    if v1 not in G.nodes() and v2 not in G.nodes():
                        G.add_edge(v1,v2,score=float(attrs[4]),linktype='tenx')
                        contigs.add(v1.split(':')[0])
                        contigs.add(v2.split(':')[0])
                        tenx_links += 1

    else:
        if args.tenx != 'abc':
            for u,v in tenx_graph.edges():
                contig_1 = u.split(':')[0]
                contig_2 = v.split(':')[0]
                if contig_1 in contig2scaffold and contig_2 in contig2scaffold:
                    scaffold_1 = contig2scaffold[contig_1]
                    scaffold_2 = contig2scaffold[contig_2]
                    test_edge(scaffold_1,scaffold_2,tenx_graph,'tenx')


'''
This function loads unitig links
'''
def load_unitig_links():
    if iteration == 1:
        if args.unitigs != 'abc':
            for u,v in unitig_graph.edges():
                if u not in G.nodes() and v not in G.nodes():
                    G.add_edge(v1,v2,score=2,linktype='unitig')
                    c1 = u.split(':')[0]
                    c2 = v.split(':')[0]
                    tiling_edges += 1
                    contigs.add(c1)
                    contigs.add(c2)

    else:
        if args.unitigs != 'abc':
            for u,v in unitig_graph.edges():
                contig_1 = u.split(':')[0]
                contig_2 = v.split(':')[0]
                if contig_1 in contig2scaffold and contig_2 in contig2scaffold:
                    scaffold_1 = contig2scaffold[contig_1]
                    scaffold_2 = contig2scaffold[contig_2]
                    test_edge(scaffold_1,scaffold_2,unitig_graph,'unitig') 


def generate_scaffold_graph():
    with open(args.directory+'/contig_links_scaled_sorted_iteration_'+str(iteration),'r') as f:
        for line in f:
            # i += 1
            # print i
            attrs = line.split()
            #print line
            if len(attrs) < 5:
                break
            v1 = attrs[0]
            v2 = attrs[1]
            c1 = attrs[0].split(':')[0]
            o1 = attrs[0].split(':')[1]
            c2 = attrs[1].split(':')[0]
            o2 = attrs[1].split(':')[1]
            #print float(attrs[4])
            all_G.add_edge(v1,v2,score=float(attrs[4]))
            #filter out links
            if iteration == 1:
                if contig_length[c1] <= int(args.cutoff) or contig_length[c2] <= int(args.cutoff):
                    continue
            if float(attrs[4]) < 1:
                # if args.unitigs !='abc':

                #     '''
                #     Before anything, check if any of the node is present in the graph. If not then no point of doing this
                #     '''
                if iteration == 1:
                    break
                    # to_continue = True
                    # for link in ["B:E", "B:B", "E:B", "E:E"]:
                    #     v1 = c1 + ':'+link[0]
                    #     v2 = c2 + ':' + link[2]
                    #     v3 = c1 + ':' + link[2]
                    #     v4 = c2 + ':' + link[0]
                    #     nodes = G.nodes()
                    #     if v1 in nodes or v2 in nodes or v3 in nodes or v4 in nodes:
                    #         to_continue = False

                    # if not to_continue:
                    #     #print 'continue'
                    #     continue


            #         if int(args.iteration) == 1:
                     

            #             '''
            #             First check if the edge is present in the unitig graph in both the orientations
            #             If the unitig graph can not provide this information, then check the GFA graph
            #             '''
            #             added = False
            #             for link in ["B:E", "B:B", "E:B", "E:E"]:
            #                 edgeType = link.split(':')
            #                 if unitig_graph.has_edge(c1+':'+edgeType[0], c2+':'+edgeType[1]):
            #                     v1 = c1+':'+edgeType[0]
            #                     v2 = c2+':'+edgeType[1]
            #                     if v1 not in G.nodes() and v2 not in G.nodes():
            #                         G.add_edge(v1,v2,score=2)
            #                         tiling_edges += 1
            #                         contigs.add(c1)
            #                         contigs.add(c2)
            #                         added = True

                                 
            #                 if unitig_graph.has_edge(c1+':'+edgeType[1], c2+':'+edgeType[0]):
            #                     v1 = c1+':'+edgeType[1]
            #                     v2 = c2+':'+edgeType[0]
            #                     if v1 not in G.nodes() and v2 not in G.nodes():
            #                         G.add_edge(v1,v2,score=2)
            #                         tiling_edges += 1
            #                         contigs.add(c1)
            #                         contigs.add(c2)
            #                         added = True
            #             if added:
            #                 print 'added'
            #             # else:
            #             #     print 'no'

            #             # if not added:
            #             #     minpath,ori,ratio = get_best_path(c1,c2)
            #             #     if minpath != -1 and ratio > 0:
            #             #         #if v1 not in G.nodes() and v2 not in G.nodes():
            #             #         print 'adding'
            #             #         gfa_edges += 1
            #             #         G.add_edge(c1+':'+ori[0],c1+':'+ori[1],score=1.5)
            #             #         contigs.add(c1)
            #             #         contigs.add(c2)

            #         #             if edgeType[0] == 'B':
            #         #                 or1 = 'REV'
            #         #             else:
            #         #                 or1 = 'FOW'

            #         #             if edgeType[1] == 'B':
            #         #                 or2 = 'FOW'
            #         #             else:
            #         #                 or2 = 'REV'
            #         #             start = c1+'$'+or1
            #         #             end = c2 + '$'+or2
            #         #             p = -1
            #         #             if start in shortest_paths:
            #         #                 if end in shortest_paths[start]:
            #         #                     p = shortest_paths[start][end]
            #         #             if minpath == -1 or p < minpath:
            #         #                 if minpath != -1:
            #         #                     ratio = p*1.0/minpath
            #         #                 minpath = p
            #         #                 ori = link
            #         #         #print ratio
            #         #         if minpath != -1 and ratio > 0:
            #         #             if v1 not in G.nodes() and v2 not in G.nodes():
            #         #                 #print 'adding'
            #         #                 gfa_edges += 1
            #         #                 G.add_edge(c1+':'+ori[0],c1+':'+ori[1],score=1.5)
            #         #                 contigs.add(c1)
        #         #                 contigs.add(c2)
                # else:
                #     '''
                #     Here too, first check 10x links and then unitigs
                #     '''

                #     if args.unitigs != 'abc':
                #         '''
                #         This is tricky. Scaffold is the series of contigs. Check just for the terminal contigs. Look up for their
                #         path in the unitig tiling and decide if to put an edge or not
                #         '''
                #         scaffold_first = previous_scaffolds[c1]
                #         scaffold_second = previous_scaffolds[c2]
                #         #if len(scaffold_first) > 2 and len(scaffold_second) > 2:
                #         first_first = scaffold_first[0]
                #         last_first = scaffold_first[-1]
                #         first_second = scaffold_second[0]
                #         last_second = scaffold_second[-1]
                        
                #         if unitig_graph.has_edge(last_first,first_second):
                #             v1 = c1+':E'
                #             v2 = c2 + ':B'
                #             if v1 not in G.nodes() and v2 not in G.nodes():
                #                 G.add_edge(v1,v2,score=2)

                #         if unitig_graph.has_edge(reverse_complement(first_first),first_second):
                #             v1 = c1 + ':B'
                #             v2 = c2 + ':B'
                #             if v1 not in G.nodes() and v2 not in G.nodes():
                #                 G.add_edge(v1,v2,score=2)

                #         if unitig_graph.has_edge(last_first,reverse_complement(last_second)):
                #             v1 = c1 + ':E'
                #             v2 = c2 + ':E'
                #             if v1 not in G.nodes() and v2 not in G.nodes():
                #                 G.add_edge(v1,v2,score=2)

                #         if unitig_graph.has_edge(reverse_complement(first_first),reverse_complement(last_second)):
                #             v1 = c1 + ':E'
                #             v2 = c2 + ':B'
                #             if v1 not in G.nodes() and v2 not in G.nodes():
                #                 G.add_edge(v1,v2,score=2)
                #                 contigs.add(c1)
                #                 contigs.add(c2)
                #                 tiling_edges += 1

                
            #     continue

            # else:
            if v1 not in G.nodes() and v2 not in G.nodes():
                #print line
                hic_edges += 1
                G.add_edge(v1,v2,score=float(attrs[4]),linktype='hic')
            # else:
            #     print "UNUSED "+ line
                contigs.add(c1)
                contigs.add(c2)
    print >> sys.stderr,  'Finished loading Hi-C links, Loading unitig links now..'


         




    '''
    Now try to add 10x and unitig links to the graph. Note that current preference is
    10x links first and then unitig tiling. 
    TODO: probably give an option to provide preference
    '''

    load_tenx_links()
    load_unitig_links()

    #Now do usual layout
    for ctg in list(contigs):
        G.add_edge(ctg+":B", ctg+":E", t="c", score=0, linktype='implicit')

    print >> sys.stderr, 'Hybrid scaffold graph loaded, nodes = ' + str(len(G.nodes())) + ' edges = ' + str(len(G.edges())) 
    print >> sys.stderr, 'Hi-C implied edges = ' + str(hic_edges)
    print >> sys.stderr, 'Unitig tiling implied edges = ' + str(tiling_edges)
    print >> sys.stderr, 'Assembly graph implied edges = ' + str(gfa_edges)
    print >> sys.stderr, '10x implied edges = ' + str(tenx_links)

'''
This function generates seeds scaffolds from the hybrid Graph G
'''
def get_seed_scaffold():

    g_idx = 1
    seed_scaffolds = {} #this stores initial long scaffolds
    to_merge = set()

    for subg in nx.connected_component_subgraphs(G):
        p = []
        for node in subg.nodes():
            if subg.degree(node) == 1:
                p.append(node)

        #If this is 2 then we have found the path!
        if len(p) == 2:
            path = nx.shortest_path(subg,p[0],p[1])
            seed_scaffolds[g_idx] = path
            g_idx += 1


        #else try to insert these contigs in the long scaffolds generated previously
        else:
            for node in subg.nodes():
                to_merge.add(node.split(':')[0])

    return seed_scaffolds, to_merge



'''
Given a small contig, it finds best scaffold where small contig can go in 
'''

def assign_small_to_seed(to_merge, seed_scaffolds):
    assignment = {}

#print to_merge

    for contig in to_merge:
        max_sum = -1
        max_path = -1
        for key in seed_scaffolds:
            path = seed_scaffolds[key]
            cur_sum = 0
            cnt = 0
            five_prime = contig+':B'
            three_prime = contig+':E'
            for node in path:
                if all_G.has_edge(five_prime,node):
                    cur_sum += all_G[five_prime][node]['score']
                    cnt += 1
                if all_G.has_edge(three_prime,node):
                    cur_sum += all_G[three_prime][node]['score']
                    cnt += 1
            if cnt != 0 and cur_sum > max_sum:
                max_sum = cur_sum
                max_path = key

        if max_sum != -1:
            assignment[contig] = max_path

    return assignment


'''
Given assignment of small to seed, this methods tries to put small scaffolds on seed in 
all possible orientation and orderings
'''
def insert(assignment, seed_scaffolds):
    to_add_later = set()
    for contig in assignment:
        #print contig
        #print assignment[contig]
        path = seed_scaffolds[assignment[contig]]
        #print path
        five_prime = contig+':B'
        three_prime = contig+':E'
        total_max = -1
        orientation = ''
        pos = -1
        #first check at all middle positions
        for i in xrange(1,len(path)-1,2):
            score_fow = -1
            score_rev = -1
            if all_G.has_edge(five_prime,path[i]) and all_G.has_edge(three_prime,path[i+1]):
                score_fow = all_G[five_prime][path[i]]['score'] + all_G[three_prime][path[i+1]]['score']
            if all_G.has_edge(three_prime,path[i]) and all_G.has_edge(five_prime,path[i+1]):
                score_rev = all_G[three_prime][path[i]]['score'] + all_G[five_prime][path[i+1]]['score']

            # print score_fow, score_rev
            if score_fow >= score_rev:
                if score_fow > total_max:
                    total_max = score_fow
                    orientation = 'fow'
                    pos = i
                else:
                    if score_rev > total_max:
                        total_max = score_rev
                        orientation = 'rev'
                        pos = i

        if total_max != -1:
            #print contig
            if orientation == 'fow':
                #print "INSERTING " + contig + ' between ' + path[pos-1] + ' and ' + path[pos+1]
                path.insert(pos+1,five_prime)
                path.insert(pos+2,three_prime)

            else:
                #print "INSERTING " + contig + ' between ' + path[pos-1] + ' and ' + path[pos+1]
                path.insert(pos+1,three_prime)
                path.insert(pos+2,five_prime)
            seed_scaffolds[assignment[contig]] = path


        else:
            #print contig
            to_add_later.add(contig)

    return seed_scaffolds

#print seed_scaffolds

#print to_merge

#assign a backbone scaffold to each small contig

# def generate_new_links(scaffolds):
#     #first generate scaffold to length
#     scaff2length = {}
#     offset = {}

#     scaffold_re = {}
#     total_scaffold_re = {}
#     total_scaffold_links = {}
#     contig_re = {}
#     with open(args.directory+'/re_counts_iteration_1','r') as f:
#         for line in f:
#             attrs = line.split()
#             contig_re[attrs[0]] = (int(attrs[1]),int(attrs[2]))

    

#     for key in scaffolds:
#         scaffold = scaffolds[key]
#         #print scaffold
#         offset[key] = {}
#         length = 0
#         try:
#             for i in xrange(1,len(scaffold),2):
#                 #   print scaffold[i]
#                 offset[key][scaffold[i].split(':')[0]] = [length]
#                 length += contig_length[scaffold[i].split(':')[0]]
#                 offset[key][scaffold[i].split(':')[0]].append(length)
#             scaff2length[key] = length
#         except:
#             print scaffold
#     #print offset
#     for key in scaffolds:
#         total_scaffold_re[key] = [0,0]
#         scaffold = scaffolds[key]
#         left = 0
#         right = 0
#         for i in xrange(1,len(scaffold),2):
#             length = scaff2length[key]
#             total_scaffold_re[key][0] += contig_re[scaffold[i].split(':')[0]][0]
#             total_scaffold_re[key][1] += contig_re[scaffold[i].split(':')[0]][1]
#             # if i < len(scaffold) - 1:
#             #     if offset[scaffold[i].split(':')[0]] <= length/3 and offset[scaffold[i+1].split(':')[0]] >= 2*length/3:
#             #         print "ambiguous"
#             if offset[key][scaffold[i].split(':')[0]][0] <= length/2 and offset[key][scaffold[i].split(':')[0]][0] <= length/2:
#                 left += contig_re[scaffold[i].split(':')[0]][0]
#             if offset[key][scaffold[i].split(':')[0]][0] >= length/2 and offset[key][scaffold[i].split(':')[0]][0] >= length/2:
#                 #print 'here'
#                 right += contig_re[scaffold[i].split(':')[0]][1]
        
#         scaffold_re[key] = [left,right]
#     #print scaffold_re
    
#     counts = nx.Graph() 
#     with open(args.directory+'/contig_links_raw_iteration_1','r') as f:
#         for line in f:
#             attrs = line.split()
#             counts.add_edge(attrs[0],attrs[1],weight=int(attrs[2]))
#     new_links = {}
#     new_links_re = {}
#     #TODO LOAD LINKS FROM THE GRAPH. ASSUME FOR NOW THAT YOU HVE LINKS
#     empty = 0
#     for key1 in scaffolds:
#         for key2 in scaffolds:
#             if key1 < key2:
#                 scaffold_1 = scaffolds[key1]
#                 scaffold_2 = scaffolds[key2]
#                 scaffold_1_length = len(scaffold_1)
#                 scaffold_2_length = len(scaffold_2)
#                 #calculate links for all possible orientations
#                 for i in xrange(1,len(scaffold_1),2):
#                     for j in xrange(1,len(scaffold_2),2):

#                         contig_1 = scaffold_1[i].split(':')[0]
#                         contig_2 = scaffold_2[j].split(':')[0]
#                         if counts.has_edge(contig_1,contig_2):
#                             new_key = key1+'$'+key2
#                             if new_key not in total_scaffold_links:
#                                 total_scaffold_links[new_key] = 0
#                             total_scaffold_links[new_key] = counts[contig_1][contig_2]['weight']
#                         if contig_1 in offset[key1] and contig_2 in offset[key2]:
#                             offset_1 = offset[key1][contig_1]
#                             offset_2 = offset[key2][contig_2]
#                             link_type = ''
#                             if offset_1[1]<= scaffold_1_length and offset_2[1] <= scaffold_2_length/3:
#                                 link_type = 'BB'
#                             if offset_1[1]<= scaffold_1_length and offset_2[0] >= 2*scaffold_2_length/3:
#                                 link_type = "BE"
#                             if offset_1[0] >= 2*scaffold_1_length/3 and offset_2[0] >= 2*scaffold_2_length/3:
#                                 link_type = "EE"
#                             if offset_1[0] >= 2*scaffold_1_length/3 and offset_2[1] <= scaffold_2_length/3:
#                                 link_type = "EB"

#                             if link_type == '':
#                                 empty += 1
#                             if counts.has_edge(contig_1,contig_2) and link_type != '':
#                                 link_key = key1+':'+link_type[0] + '$'+key2+':'+link_type[1]
#                                 if link_key not in new_links:
#                                     new_links[link_key] = 0
#                                     new_links_re[link_key] = 0
#                                 new_links[link_key] += counts[contig_1][contig_2]['weight']
#                                 if link_type == 'BB':
#                                     new_links_re[link_key] += scaffold_re[key1][0] + scaffold_re[key2][0]
#                                 if link_type == 'BE':
#                                     new_links_re[link_key] += scaffold_re[key1][0] + scaffold_re[key2][1]
#                                 if link_type == 'EB':
#                                     new_links_re[link_key] += scaffold_re[key1][1] + scaffold_re[key2][0]
#                                 if link_type == 'EE':
#                                     new_links_re[link_key] += scaffold_re[key1][1] + scaffold_re[key2][1]

#     print "Confused links = "+str(empty)
#     ofile = open(args.directory+'/contig_links_iteration_'+str(iteration+1),'w')
#     for link in new_links:
#         contigs = link.split('$')
#         if new_links_re[link] == 0:
#             new_links_re[link] = 1
#         ofile.write(contigs[0]+'\t'+contigs[1]+'\t'+ str(new_links[link]*1.0/new_links_re[link])+'\t'+str(new_links[link])+'\n')

#     ofile.close()
#     ofile = open(args.directory+'/re_counts_iteration_'+str(iteration+1),'w')
#     for key in total_scaffold_re:
#         ofile.write(key+'\t'+str(total_scaffold_re[key][0])+'\t'+str(total_scaffold_re[key][1])+'\n')
#     ofile.close()

#     ofile = open(args.directory+'/contig_links_raw_iteration_'+str(iteration+1),'w')
#     for key in total_scaffold_links:
#         ofile.write(key.split('$')[0]+'\t'+key.split('$')[1]+'\t'+str(total_scaffold_links[key])+'\n')

def update_bed(expanded_scaffold):
    contig2scaffold = {}
    contig2info = {}
    scaffold_length = {}

    re_counts = {}
    with open(args.directory+'/re_counts_iteration_1','r') as f:
        for line in f:
            attrs = line.split()
            re_counts[attrs[0]] = (int(attrs[1]),int(attrs[2]))

    #print re_counts

    for key in expanded_scaffold:
        path = expanded_scaffold[key]
        scaffold_length[key] = 0
        offset = 0
        for i in xrange(0,len(path),2):
            contig = path[i].split(':')[0]
            contig2scaffold[contig] = key
            ori = path[i].split(':')[1] + path[i+1].split(':')[1]
            if ori == 'BE':
                contig2info[contig] = (offset,'FOW')
            else:
                contig2info[contig] = (offset,'REV')
            offset += contig_length[contig]
            scaffold_length[key] += contig_length[contig]

    scaffold_re = {}
    for key in expanded_scaffold:
        path = expanded_scaffold[key]
        length = scaffold_length[key]
        offset = 0
        s_left = 0
        s_right = 0
        if len(path) == 2:
            contig = path[0].split(':')[0]
            scaffold_re[key] = re_counts[contig]
        else:
            for i in xrange(0,len(path),2):
                contig = path[i].split(':')[0]
                contig2scaffold[contig] = key
                left,right = re_counts[contig]
                if offset <= length/2 and i+2 < len(path):
                    if contig2info[path[i+2].split(':')[0]][0] <= length/2:
                        s_left += (left + right) 
                    else:                
                        contig_next = path[i+2].split(':')[0]
                        if contig2info[contig_next][0] >= length/2:
                            left_part = length/2 - offset
                            right_part = contig2info[path[i+2].split(':')[0]][0] - length/2
                            s_left += left*left_part/contig_length[contig]
                            s_right += right*right_part/contig_length[contig]

                if offset <= length/2 and i + 2 >= len(path):
                    left_part = length/2 - offset
                    right_part = length/2
                    s_left += left*left_part/contig_length[contig]
                    s_right += right*right_part/contig_length[contig]

                if offset >= length/2:
                    s_right += (left+right)
                offset += contig_length[contig]
                #scaffold_length[key] += contig_length[contig]
            scaffold_re[key] = (s_left,s_right)
    
    o_lines = ""
    count = 0
    if not os.path.isfile(args.directory+'/alignment_iteration_'+str(iteration+1)+'.bed'):
        output = open(args.directory+'/alignment_iteration_'+str(iteration+1)+'.bed','w')
        with open(args.directory+'/alignment_iteration_1.bed','r') as f:
            for line in f:
                attrs = line.split()
                contig = attrs[0]
                start = int(attrs[1])
                end = int(attrs[2])
                if contig in contig2scaffold:
                    scaffold = contig2scaffold[contig]
                    info = contig2info[contig]
                    new_start = start + info[0]
                    new_end = end + info[0]
                    o_lines += scaffold+'\t'+str(new_start)+'\t'+str(new_end)+'\t'+attrs[3]+'\n'
                    count += 1
                if count == 100000:
                    output.write(o_lines)
                    count = 0
                    o_lines = ""

            #write remaining lines
            output.write(o_lines)
            output.close()
    len_output = open(args.directory+'/scaffold_length_iteration_'+str(iteration+1),'w')
    for key in scaffold_length:
        len_output.write(key+'\t'+str(scaffold_length[key])+'\n')
    len_output.close()

    
    re_out = open(args.directory+'/re_counts_iteration_'+str(iteration+1),'w')
    for key in scaffold_re:
        re_out.write(key+'\t'+str(scaffold_re[key][0])+'\t'+str(scaffold_re[key][1])+'\n')
    re_out.close()




#print assignment
#now try to place these contigs on these paths in all possible orientations and at all possible positions
#print assignment


def merge(contigs):
    #print "MERGING " + str(contigs)
    scaffolds = []
    subg = all_G.subgraph(contigs)
    best_hic_graph = nx.Graph()
    edges = []
    contigs = set()
    for u,v,data in subg.edges(data=True):
        edges.append((u,v,data['score']))
    edges.sort(key=lambda x: x[2],reverse=True)
    #print edges
    for u,v,score in edges:
        if u not in best_hic_graph.nodes() and v not in best_hic_graph.nodes():
            best_hic_graph.add_edge(u,v,score=score)
            contigs.add(u.split(':')[0])
            contigs.add(v.split(':')[0])
    
    for contig in contigs:
        best_hic_graph.add_edge(contig+':B',contig+':E',score=0)

    
    #print best_hic_graph.nodes()

    for g in nx.connected_component_subgraphs(best_hic_graph):
        #print g.nodes()
        p = []
        for node in g.nodes():
            if g.degree(node) == 1:
                p.append(node)

        if len(p) == 2:
            path = nx.shortest_path(g,p[0],p[1])
            scaffolds.append(path)
        else:
            contigs = set()
            for node in g.nodes():
                contigs.add(node.split(':')[0])
            for each in contigs:
                scaffolds.append([each+':B',each+':E'])
    #print scaffolds
    return scaffolds


'''
Call all the methods here
'''

generate_scaffold_graph()
seed_scaffolds,to_merge =  get_seed_scaffold()
assignment = assign_small_to_seed(to_merge,seed_scaffolds)
seed_scaffolds = insert(assignment,seed_scaffolds)

#try to merge alternating scaffolds
merged = {}
final_scaffolds = {}
scaffold_id = 1
# scaffold_id = len(previous_scaffolds) + 1 #important to have this adjustment to keep all scaffold names unique
# for key1 in seed_scaffolds:
#     for key2 in seed_scaffolds:
#         if key1 not in merged and key2 not in merged:
#             if key1 < key2:
#                 count = 0
#                 path1 = seed_scaffolds[key1]
#                 path2 = seed_scaffolds[key2]
#                 for i in xrange(len(path1)):
#                     for j in xrange(len(path2)):
#                         if all_G.has_edge(path1[i],path2[j]):
#                             if all_G[path1[i]][path2[j]]['score'] >= 0.1:
#                                 count += 1

#                 percent = count*100.0/(len(path1)*len(path2))
#                 #print "Percent = "+ str(percent)

#                 if percent >= 80:
#                     merged[key1] = True
#                     merged[key2] = True
#                     contigs = set()
#                     for each in path1:
#                         contigs.add(each)
#                     for each in path2:
#                         contigs.add(each)

#                     merged_scaffolds = merge(contigs)

#                     for each in merged_scaffolds:
#                         #print each
#                         final_scaffolds[scaffold_id] = each
#                         scaffold_id += 1

#add scaffolds that are not merged
for key in seed_scaffolds:
    if key not in merged:
        final_scaffolds[scaffold_id] = seed_scaffolds[key]
        scaffold_id += 1

#print to_add_later
for key in to_add_later:
    final_scaffolds[scaffold_id] = [key+':B',key+':E']
    scaffold_id += 1 
#now here expand the scaffolds with the map loaded before. Assumption here is that the map will not contain any path that has scaffolds in it.
#Everything is contig

scaffolded_contigs = {}
expanded_scaffold_paths = {}
visited = {}
for key in final_scaffolds:
    path = final_scaffolds[key]
    new_path = []
    #print path
    for i in xrange(0,len(path)-1,2):
        #print path
        contig = path[i].split(':')[0]
        scaffolded_contigs[contig] = True
        if contig[0] == 's':
            visited[contig] = True
            #print contig
            actual_path = previous_scaffolds[contig]
            link_type = path[i].split(':')[1] + path[i+1].split(':')[1]
            if link_type != 'BE':
              actual_path = actual_path[::-1]
            #print actual_path
            for each in actual_path:
              new_path.append(each)
        else:
          new_path.append(path[i])
          new_path.append(path[i+1])


    expanded_scaffold_paths['scaffold_'+str(key)] = new_path


#check scaffolds which were not scaffolded in this round but present in previous round

if iteration > 1:
    for key in previous_scaffolds:
        #key1 = 'scaffold_'+str(key)
        if key not in scaffolded_contigs:
            #print key
            expanded_scaffold_paths['scaffold_'+str(scaffold_id)] = previous_scaffolds[key]
            scaffold_id += 1 


if iteration == 1:
    for contig in contig_length:
        if contig not in scaffolded_contigs:
            #print contig
            expanded_scaffold_paths["scaffold_"+str(scaffold_id)] = [contig+':B',contig+':E']
            scaffold_id += 1


new_id = 1
new_scaffolds = {}
for key in expanded_scaffold_paths:
    new_scaffolds['scaffold_'+str(new_id)] = expanded_scaffold_paths[key]
    new_id += 1

expanded_scaffold_paths = new_scaffolds
new_scaffolds = {}      

pickle.dump(expanded_scaffold_paths,open(args.directory+'/scaffolds_iteration_'+str(args.iteration)+'.p','w'))
update_bed(expanded_scaffold_paths)



