import os
import argparse
import sys
import subprocess
from subprocess import Popen, PIPE

def check(path):
	with open(path,'r') as f:
		for line in f:
			attrs = line.split()
			if float(attrs[4]) >= 1:
				return False
			else:
				return True

def main():
	#bin=os.path.dirname(os.path.abspath(__file__))
	parser = argparse.ArgumentParser(description="SALSA Iterative Pipeline")
	parser.add_argument('-a','--assembly',help='Path to initial assembly',required=True)
	parser.add_argument('-l','--length',help='Length of contigs at start',required=True)
	parser.add_argument('-b','--bed',help='Bed file of alignments sorted by read names',required=True)
	parser.add_argument('-o','--output',help='Output directory to put results',required=False,default='SALSA_output')
	parser.add_argument('-c','--cutoff',help='Minimum contig length to scaffold',required=False,default=1000)
	parser.add_argument('-g','--gfa',help='GFA file for assembly',required=False,default='abc')
	parser.add_argument('-u','--unitigs',help='The tiling of unitigs to contigs in bed format',required=False,default='abc')
	parser.add_argument('-t','--tenx',help='10x links tab separated file',required=False,default='abc')
	parser.add_argument('-e','--enzyme',help='Restriction Enzyme used for experiment',required=False,default='AAGCTT')


	args = parser.parse_args()


	#iteration counter
	iter_num = 1

	if not os.path.exists(args.output):
		os.mkdir(args.output)

	os.system('cp '+args.length+' '+args.output+'/scaffold_length_iteration_1')
	if not os.path.isfile(args.output+'/alignment_iteration_1.bed'):
		#os.system('cp '+args.bed+' '+args.output+'/alignment_iteration_1.bed')
		os.symlink(args.bed,args.output+'/alignment_iteration_1.bed')
	log = open(args.output+'/commands.log','w')

	#First get RE sites
	if not os.path.isfile(args.output+'/re_counts_iteration_'+str(iter_num)):
		try:
			cmd = 'python RE_sites.py -a '+args.assembly + ' -e '+ args.enzyme + ' > '+ args.output+'/re_counts_iteration_'+str(iter_num)
			print cmd
			log.write(cmd+'\n')
			p = subprocess.check_output(cmd,shell=True)

		except subprocess.CalledProcessError as err:
			# if os.path.isfile(args.output+'/re_counts_iteration_'+str(iter_num)):
			# 	os.system(args.output+'/RE_sites_iteration_'+str(iter_num))
			print >> sys.stderr, str(err.output)
			sys.exit(1)

	#Now compute normal links with old new_links code
	print >> sys.stderr, "Starting Iteration "+ str(iter_num)
	if not os.path.isfile(args.output+'/contig_links_iteration_'+str(iter_num)):
		try:
			cmd = 'python new_links_test.py -b '+ args.bed + ' -d '+ args.output +' -i '+str(iter_num)
			print cmd
			log.write(cmd+'\n')
			p = subprocess.check_output(cmd,shell=True)

		except subprocess.CalledProcessError as err:
			print >> sys.stderr, str(err.output)
			if os.path.isfile(args.output+'/contig_links_iteration_'+str(iter_num)):
				os.system('rm '+args.output+'/contig_links_iteration_'+str(iter_num))
			sys.exit(1)

	#now use Serge's code to calculate 
	if not os.path.isfile(args.output+'/contig_links_scaled_iteration_'+str(iter_num)):
		try:
			cmd =  'python fast_scaled_scores.py -d '+args.output+' -i '+str(iter_num)
			log.write(cmd+'\n')
			print cmd
			p = subprocess.check_output(cmd,shell=True)
		except subprocess.CalledProcessError as err:
			if os.path.isfile(args.output+'/contig_links_scaled_iteration_'+str(iter_num)):
				os.system('rm '+args.output+'/contig_links_scaled_iteration_'+str(iter_num))
			print >> sys.stderr, str(err.output)
			sys.exit(1)

	#Sort the links by column 5 
	if not os.path.isfile(args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num)):
		try:
			cmd = 'sort -k 5 -gr '+args.output+'/contig_links_scaled_iteration_'+str(iter_num)+ ' > '+ args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num)
			log.write(cmd+'\n')
			print cmd
			p = subprocess.check_output(cmd,shell=True)
		except subprocess.CalledProcessError as err:
			if os.path.isfile(args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num)):
				os.system('rm '+args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num))
			print >> sys.stderr, str(err.output)
			sys.exit(1)

	if not os.path.isfile(args.output+'/scaffolds_iteration_1.p'):
		try:
			cmd = 'python layout_refactory_unitigs_test.py -x '+args.gfa + ' -u '+args.unitigs+' -t '+args.tenx+' -l '+args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num) +' -c ' +args.cutoff+' -i 1 -d '+args.output
			log.write(cmd+'\n')
			print cmd
			p = subprocess.check_output(cmd,shell=True)

		except subprocess.CalledProcessError as err:
			if os.path.isfile(args.output+'scaffolds_iteration_1.p'):
				os.system('rm '+args.output+'scaffolds_iteration_1.p')
			print >> sys.stderr, str(err.output)
			sys.exit(1)

	iter_num += 1

	#now do iterative
	while True:
		print >> sys.stderr, "Starting Iteration "+ str(iter_num)
		if not os.path.isfile(args.output+'/contig_links_iteration_'+str(iter_num)):
			try:
				cmd = 'python new_links_test.py -b '+ args.output+'/alignment_iteration_'+str(iter_num)+'.bed' + ' -d '+ args.output +' -i '+str(iter_num)
				log.write(cmd+'\n')
				p = subprocess.check_output(cmd,shell=True)
				os.system("rm "+args.output+'/alignment_iteration_'+str(iter_num)+'.bed')

			except subprocess.CalledProcessError as err:
				print >> sys.stderr, str(err.output)
				if os.path.isfile(args.output+'/contig_links_iteration_'+str(iter_num)):
					os.system('rm '+args.output+'/contig_links_iteration_'+str(iter_num))
				sys.exit(1)

		print >> sys.stderr, "Starting Iteration "+ str(iter_num)
		if not os.path.isfile(args.output+'/contig_links_scaled_iteration_'+str(iter_num)):
			try:
				cmd =  'python fast_scaled_scores.py -d '+args.output+' -i '+str(iter_num)
				log.write(cmd+'\n')
				p = subprocess.check_output(cmd,shell=True)
			except subprocess.CalledProcessError as err:
				if os.path.isfile(args.output+'/contig_links_scaled_iteration_'+str(iter_num)):
					os.system('rm '+args.output+'/contig_links_scaled_iteration_'+str(iter_num))
				print >> sys.stderr, str(err.output)
				sys.exit(1)

		if not os.path.isfile(args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num)):
			try:
				cmd = 'sort -k 5 -gr '+args.output+'/contig_links_scaled_iteration_'+str(iter_num)+ ' > '+ args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num)
				log.write(cmd+'\n')
				p = subprocess.check_output(cmd,shell=True)
			except subprocess.CalledProcessError as err:
				if os.path.isfile(args.output+'/new_links_scaled_sorted_iteration_'+str(iter_num)):
					os.system('rm '+args.output+'/new_links_scaled_sorted_iteration_'+str(iter_num))
				print >> sys.stderr, str(err.output)
				sys.exit(1)
		#NOW check if any useful link here
		if check(args.output+'/contig_links_scaled_sorted_iteration_'+str(iter_num)):
			break

		if not os.path.isfile(args.output+'scaffolds_iteration_'+str(iter_num)+'.p'):
			try:
				cmd = 'python layout_refactory_unitigs_test.py -x abc -u '+args.unitigs+' -t '+args.tenx+' -l '+args.output+'/new_links_scaled_sorted_iteration_'+str(iter_num) +' -c ' +args.cutoff+' -i '+str(iter_num)+' -d '+args.output
				log.write(cmd+'\n')
				p = subprocess.check_output(cmd,shell=True)

			except subprocess.CalledProcessError as err:
				if os.path.isfile(args.output+'scaffolds_iteration_'+str(iter_num)+'.p'):
					os.system('rm '+args.output+'scaffolds_iteration_'+str(iter_num)+'.p')
				print >> sys.stderr, str(err.output)
				sys.exit(1)


		

		iter_num += 1

	log.close()

if __name__ == '__main__':
	main()

