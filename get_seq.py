from pbcore.io import FastaReader
from pbcore.io import FastaWriter
import pickle
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-a","--cleaned",help="cleaned assembly")
parser.add_argument("-f","--scaffold",help="final scaffold file")
parser.add_argument("-g","--agp",help="agp file")
parser.add_argument("-p","--map",help="pickle map of scaffolds")

args = parser.parse_args()

revcompl = lambda x: ''.join([{'A':'T','C':'G','G':'C','T':'A','N':'N','R':'N','M':'N','Y':'N','S':'N','W':'N','K':'N','a':'t','c':'g','g':'c','t':'a','n':'n',' ':'',}[B] for B in x][::-1])
scaff_map = pickle.load(open(args.map,'r'))

contig_length = {}
id2seq = {}

f = FastaReader(args.cleaned)
for record in f:
    id2seq[record.id] = record.sequence
    contig_length[record.id] = len(record.sequence)


c_id = 1
line = ""
agp_output = open(args.agp,'w')
writer = FastaWriter(args.scaffold)
for key in scaff_map:
    start = 1
    local_comp = 1
    if len(scaff_map[key]) >= 4:
        path = scaff_map[key]
        scaff_len = 0
        curr_contig = ""
        #print c_id
        for i in range(0,len(path)-1,2):
            line += 'scaffold_'+str(c_id)
            line += '\t'
            line += str(start)
            line += str('\t')
            curr = path[i]
            next = path[i+1]
            curr = curr.split(':')
            next = next.split(':')
            curr_len = contig_length[curr[0]]
            scaff_len += curr_len
            end = curr_len + start - 1
            line += str(end)
            line += '\t'
            start = end + 1
            line += str(local_comp)
            local_comp += 1
            line += ('\tW\t' + curr[0] + '\t' + '1\t')
            line += str(curr_len)
            line += '\t'
            #print curr
            if curr[1] == 'B' and next[1] == 'E':
                curr_contig += id2seq[curr[0]]
                line += '+\t'
            if curr[1] == 'E' and next[1] == 'B':
                line += '-\t'
                #print id2seq[curr[0]]
                curr_contig += revcompl(id2seq[curr[0]])

            agp_output.write(line+'\n')
            line = ""
            line += 'scaffold_'+str(c_id)+'\t'
            line += str(start) +'\t'
            end = 500 + start - 1
            line += str(end) + '\t'
            start = end + 1
            line += str(local_comp) +'\t'
            local_comp += 1
            line += 'N\t500\tscaffold\tyes\tna'
            agp_output.write(line+'\n')
            line="" 
            if i != len(path) - 2:
                for j in range(0,500):
                    curr_contig += 'N'
        # rec = SeqRecord(Seq(curr_contig,generic_dna),id='scaffold_'+str(c_id))
        # recs.append(rec)
        #print c_id
        writer.writeRecord('scaffold_'+str(c_id),curr_contig)
        c_id += 1