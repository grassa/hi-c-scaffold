import argparse
import pickle

parser = argparse.ArgumentParser()
parser.add_argument('-p','--map', help='Pickle map of assembly in previous stage', required=False)
parser.add_argument('-l','--length', help='Length of contigs at start', required=True)
parser.add_argument('-i','--iteration', help='Iteration number', required=False)
parser.add_argument('-r','--resites', help='cut sites at start', required=False)
parser.add_argument('-o','--olen', help='length output', required=False)
parser.add_argument('-x','--reoutput', help='resites output', required=False)
args = parser.parse_args()

contig_length = {}
cut_sites = {}

with open(args.length,'r') as f:
	for line in f:
		attrs = line.split()
		contig_length[attrs[0]] = int(attrs[1])

with open(args.resites,'r') as f:
	for line in f:
		attrs = line.split()
		cut_sites[attrs[0]] = [int(attrs[1]),int(attrs[2])]

scaffold_map = pickle.load(open(args.map,'r'))

lenfile = open(args.olen,'w')
resitefile = open(args.reoutput,'w')
for scaffold in scaffold_map:
	key = scaffold
	length = 0
	cut_left = 0
	cut_right = 0
	path = scaffold_map[scaffold]
	for i in xrange(0,len(path)-1,2):
		contig = path[i].split(':')[0]
		length += contig_length[contig]
		cut_left += cut_sites[contig][0]
		cut_right += cut_sites[contig][1]

	lenfile.write(key+'\t'+str(length)+'\n')
	resitefile.write(key+'\t'+str(cut_left)+'\t'+str(cut_right)+'\n')


