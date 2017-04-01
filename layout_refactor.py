import networkx as nx
import sys
import operator
import pickle
import argparse

parser = argparse.ArgumentParser()
#parser.add_argument('-a','--assembly', help='Contig assembly', required=False)
parser.add_argument('-x','--graph', help='GFA file for assembly graph', required=False,default='abc')
parser.add_argument('-l','--links', help='Links sorted by relative score', required=True)
parser.add_argument('-g','--alignment', help='Bed file for alignments', required=False)
parser.add_argument('-c','--cutoff', help='Minimum length contig to consider for scaffolding', required=False)
parser.add_argument('-i','--iteration',help='Iteration number',required=False,default=1)
parser.add_argument('-p','--map',help='Pickle map of scaffolds generated in previous iteration',required=False)
parser.add_argument('-n','--nextmap',help='Pickle map of scaffolds generated in this iteration',required=False)
parser.add_argument('-b','--bed', help='Output bed file for next iteration', required=False)
parser.add_argument('-s','--size', help='Size of previous scaffolds/contigs', required=False)
parser.add_argument('-f','--first', help='Original contig lengths', required=False)
args = parser.parse_args()


def parse_fasta(fh):
    fa = {}
    current_short_name = None
    # Part 1: compile list of lines per sequence
    for ln in fh:
        if ln[0] == '>':
            # new name line; remember current sequence's short name
            long_name = ln[1:].rstrip()
            current_short_name = long_name.split()[0]
            fa[current_short_name] = []
        else:
            # append nucleotides to current sequence
            fa[current_short_name].append(ln.rstrip())
    # Part 2: join lists into strings
    for short_name, nuc_list in fa.iteritems():
        # join this sequence's lines into one long string
        fa[short_name] = ''.join(nuc_list)
    return fa

G = nx.Graph()
contigs = set()
all_G = nx.Graph()
contig_length = {}

contig_lengths_original = {}
with open(args.first,'r') as f:
    for line in f:
        attrs = line.split()
        contig_lengths_original[attrs[0]] = int(attrs[1])

OVl_G = nx.Graph()

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

def load_GFA(path):
    f = open(path,'r')
    for line in f:
       line = line.strip().split()

       #S   tig00000001 *   LN:i:20388
       #L   tig00023281 -   tig00008904 +   21556M  cv:A:F
       if (line[0] == "S"):
          OVl_G.add_node(line[1]+'_pilon'+"$REV", name=line[1], length=line[3].split(":")[2])
          OVl_G.add_node(line[1]+'_pilon'+"$FOW", name=line[1], length=line[3].split(":")[2])
       if (line[0] == "L"):
          #BUG BUG BUG in GFA inverting edges on input
          ctg1=line[3]+"_pilon"
          ctg2=line[1]+"_pilon"
          tmp=line[2]
          line[2]=line[4]
          line[4]=tmp

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
    f.close()



iteration = int(args.iteration)

if iteration == 1 or args.graph !='abc':
    print 'loading GFA'
    load_GFA(args.graph)

print 'DONE'

try:
    previous_scaffolds = pickle.load(open(args.map,'r'))
except:
    previous_scaffolds = {}

with open(args.size,'r') as f:
    for line in f:
        attrs = line.split()
        contig_length[attrs[0]] = int(attrs[1])

#method that updates the read mapping based on scaffolds constructed
def update_bed(scaffolds):
    contig2coords = {}
    # with open(sys.argv[3],'r') as f:
    #     for line in f:
    #         attrs = line.split()
    #         if attrs[0] not in contig2coords:
    #             contig2coords[attrs[0]] = []
    #         contig2coords[attrs[0]].append((int(attrs[1]),int(attrs[2]),attrs[3]))

    scaffolds2coords = {}
    contig_in_scaffolds = {}
    scaffolds2offset = {}
    ofile = open(args.bed,'w')
    for key in scaffolds:
        offset = {}
        scaffolds2coords[key] = []
        scaffolds2offset[key] = {}
        path = scaffolds[key]
        sum = 0
        for i in xrange(0,len(path),2):
            contig = path[i].split(':')[0]
            contig_in_scaffolds[contig] = key
            scaffolds2offset[key][contig] = sum
            #offset[contig] = sum 
            sum += contig_lengths_original[contig]

    with open(args.alignment,'r') as f:
        for line in f:
            attrs = line.split()
            if attrs[0] in contig_in_scaffolds:
                curr_scaf = contig_in_scaffolds[attrs[0]]
                offset  = scaffolds2offset[curr_scaf][attrs[0]]
                start = int(attrs[1]) + offset
                end = int(attrs[2]) + offset
                ofile.write(curr_scaf+'\t'+str(start)+'\t'+str(end)+'\t'+str(attrs[3])+'\t0\t0\n')
            else:
                ofile.write(line)

    ofile.close()
#Load the best score graph first, keep log of used and unused links
print 'LOADING LINKS'
with open(args.links,'r') as f:
    for line in f:
        attrs = line.split()
        if len(attrs) < 5:
            break
        v1 = attrs[0]
        v2 = attrs[1]
        c1 = attrs[0].split(':')[0]
        o1 = attrs[0].split(':')[1]
        c2 = attrs[1].split(':')[0]
        o2 = attrs[1].split(':')[1]
        all_G.add_edge(v1,v2,score=float(attrs[4]))
        #filter out links
        if contig_length[c1] <= int(args.cutoff) or contig_length[c2] <= int(args.cutoff):
            continue
        if float(attrs[4]) < 1:
            if iteration == 1:
                #print 'here'
                #p = get_path_two(c1+"-"+o1, c2+"-"+o2)
                #print p
                #print "Got path between %s %s %s of length %s best is %s %s"%(line[1], line[2], link, len(p), ori, minpath)
                [minpath, ori, ratio] = get_best_path(c1,c2)
                if (int(minpath) != -1 and ratio > 0):
                    #print 'here'
                    if v1 not in G.nodes() and v2 not in G.nodes():
                        print 'adding'
                        G.add_edge(v1,v2,score=1.5)
                        contigs.add(c1)
                        contigs.add(c2)
            
            continue

        if v1 not in G.nodes() and v2 not in G.nodes():
            G.add_edge(v1,v2,score=float(attrs[4]))
        # else:
        #     print "UNUSED "+ line
        contigs.add(c1)
        contigs.add(c2)

#Now do usual layout
for ctg in list(contigs):
    G.add_edge( ctg+":B", ctg+":E", t="c", score=0)

print len(G.nodes())
print len(G.edges())

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


#print seed_scaffolds

#print to_merge

#assign a backbone scaffold to each small contig

assignment = {}


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

#print assignment
#now try to place these contigs on these paths in all possible orientations and at all possible positions
#print assignment
for contig in assignment:
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
        if orientation == 'fow':
            print "INSERTING " + contig + ' between ' + path[pos-1] + ' and ' + path[pos+1]
            path.insert(pos+1,five_prime)
            path.insert(pos+2,three_prime)

        else:
            print "INSERTING " + contig + ' between ' + path[pos-1] + ' and ' + path[pos+1]
            path.insert(pos+1,three_prime)
            path.insert(pos+2,five_prime)


def merge(contigs):
    print "MERGING"
    scaffolds = []
    subg = all_G.subgraph(contigs)
    best_hic_graph = nx.Graph()
    edges = []
    contigs = set()
    for u,v,data in subg.edges(data=True):
        edges.append((u,v,data['score']))
    edges.sort(key=lambda x: x[2],reverse=True)
    for u,v,score in edges:
        if u not in best_hic_graph.nodes() and v not in best_hic_graph.nodes():
            best_hic_graph.add_edge(u,v,score=score)
            contigs.add(u.split(':')[0])
            contigs.add(v.split(':')[0])

    for contig in contigs:
        best_hic_graph.add_edge(contig+':B',contig+':E',score=0)

    for g in nx.connected_component_subgraphs(best_hic_graph):
        p = []
        for node in g.nodes():
            if g.degree(node) == 1:
                p.append(node)

        if len(p) == 2:
            path = nx.shortest_path(g,p[0],p[1])
            scaffolds.append(path)

    return scaffolds

#try to merge alternating scaffolds
merged = {}
final_scaffolds = {}
scaffold_id = len(previous_scaffolds) + 1 #important to have this adjustment to keep all scaffold names unique
for key1 in seed_scaffolds:
    for key2 in seed_scaffolds:
        if key1 not in merged and key2 not in merged:
            if key1 < key2:
                count = 0
                path1 = seed_scaffolds[key1]
                path2 = seed_scaffolds[key2]
                for i in xrange(len(path1)):
                    for j in xrange(len(path2)):
                        if all_G.has_edge(path1[i],path2[j]):
                            if all_G[path1[i]][path2[j]]['score'] >= 0.1:
                                count += 1

                percent = count*100.0/(len(path1)*len(path2))
                #print "Percent = "+ str(percent)

                if percent >= 60:
                    merged[key1] = True
                    merged[key2] = True
                    contigs = set()
                    for each in path1:
                        contigs.add(each)
                    for each in path2:
                        contigs.add(each)

                    merged_scaffolds = merge(contigs)
                    for each in merged_scaffolds:
                        final_scaffolds[scaffold_id] = each
                        scaffold_id += 1

#add scaffolds that are not merged
for key in seed_scaffolds:
    if key not in merged:
        final_scaffolds[scaffold_id] = seed_scaffolds[key]
        scaffold_id += 1

#now here expand the scaffolds with the map loaded before. Assumption here is that the map will not contain any path that has scaffolds in it.
#Everything is contig
expanded_scaffold_paths = {}
visited = {}
for key in final_scaffolds:
    path = final_scaffolds[key]
    new_path = []
    #print path
    for i in xrange(0,len(path)-1,2):
      contig = path[i].split(':')[0]
      if contig[0] == 's':
        visited[contig] = True
        actual_path = previous_scaffolds[contig]
        link_type = path[i].split(':')[0] + path[i].split(':')[1]
        if link_type == 'BE':
          continue
        else:
          actual_path = actual_path[::-1]

          for each in actual_path:
              new_path.append(each)
      else:
          new_path.append(path[i])
          new_path.append(path[i+1])

    expanded_scaffold_paths['scaffold_'+str(key)] = new_path


#check scaffolds which were not scaffolded in this round but present in previous round
for key in previous_scaffolds:
    #key1 = 'scaffold_'+str(key)
    if key not in visited:
        expanded_scaffold_paths[key] = previous_scaffolds[key]
#write expanded scaffolds as pickle file
pickle.dump(expanded_scaffold_paths,open(args.nextmap,'w'))
#update bed file
update_bed(expanded_scaffold_paths)


#Now add scaffolds from previous iteration not seen in current iteration

# for key in final_scaffolds:
#     print final_scaffolds[key]
#revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N','R':'N','M':'N','Y':'N','S':'N','W':'N','K':'N','a':'t','c':'g','g':'c','t':'a','n':'n',' ':'',}[B] for B in x][::-1])

# seq = FastaReader(sys.argv[3])
# id2seq = {}
# for each in seq:
#     id2seq[each.id] = each.sequence
# owrite = FastaWriter(sys.argv[4])
# scaff_id = 1
# for key in final_scaffolds:
#     path = final_scaffolds[key]
#     print final_scaffolds[key]
#     curr_contig = ''
#     for i in range(0,len(path)-1,2):
#         curr = path[i]
#         next = path[i+1]
#         curr = curr.split(':')
#         next = next.split(':')
#         if curr[1] == 'B' and next[1] == 'E':
#             curr_contig += id2seq[curr[0]]
#         else:
#             curr_contig += revcompl(id2seq[curr[0]])

#     owrite.writeRecor d('scaffold_'+str(scaff_id),curr_contig)

#print len(final_scaffolds)#
#pickle.dump(final_scaffolds,open('scaffolds_2.p','wb'))



