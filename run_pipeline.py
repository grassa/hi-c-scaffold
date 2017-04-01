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
	parser.add_argument('-c','--cutoff',help='Minimum contig length to scaffold',required=False,default=50000)
	parser.add_argument('-g','--gfa',help='GFA file for assembly',required=False)
	parser.add_argument('-e','--enzyme',help='Restriction Enzyme used for experiment',required=False,default='AAGCTT')


	args = parser.parse_args()

	#iteration counter
	iter_num = 1

	if not os.path.exists(args.output):
		os.mkdir(args.output)

	log = open(args.output+'/commands.log','w')

	#First get RE sites
	if not os.path.isfile(args.output+'/RE_sites_iteration_'+str(iter_num)):
		try:
			cmd = 'python RE_sites.py -a '+args.assembly + ' -e '+ args.enzyme + ' > '+ args.output+'/RE_sites_iteration_'+str(iter_num)
			log.write(cmd+'\n')
			p = subprocess.check_output(cmd,shell=True)

		except subprocess.CalledProcessError as err:
			print >> sys.stderr, str(err.output)
			sys.exit(1)

	#Now compute normal links with old new_links code
	if not os.path.isfile(args.output+'/new_links_iteration_'+str(iter_num)):
		try:
			cmd = 'python new_links.py -b '+ args.bed + ' -c '+ args.output+'/RE_sites_iteration_'+str(iter_num) + ' -l '+args.length + ' > '+ args.output+'/new_links_iteration_'+str(iter_num)
			log.write(cmd+'\n')
			p = subprocess.check_output(cmd,shell=True)

		except subprocess.CalledProcessError as err:
			print >> sys.stderr, str(err.output)
			sys.exit(1)

	#now use Serge's code to calculate 
	if not os.path.isfile(args.output+'/new_links_scaled_iteration_'+str(iter_num)):
		try:
			cmd =  'python scaled_scores.py '+args.output+'/new_links_iteration_'+str(iter_num)+ ' > ' + args.output+'/new_links_scaled_iteration_'+str(iter_num)
			log.write(cmd+'\n')
			p = subprocess.check_output(cmd,shell=True)
		except subprocess.CalledProcessError as err:
			print >> sys.stderr, str(err.output)
			sys.exit(1)

	#Sort the links by column 5 
	if not os.path.isfile(args.output+'/new_links_scaled_iteration_'+str(iter_num)+'_sorted'):
		try:
			cmd = 'sort -k 5 -gr '+args.output+'/new_links_scaled_iteration_'+str(iter_num)+ ' > '+ args.output+'/new_links_scaled_iteration_'+str(iter_num)+'_sorted'
			log.write(cmd+'\n')
			p = subprocess.check_output(cmd,shell=True)
		except subprocess.CalledProcessError as err:
			print >> sys.stderr, str(err.output)
			sys.exit(1)

	if not os.path.isfile(args.output+'/scaffolds_'+str(iter_num)+'.p'):
		try:
			cmd = 'python layout_refactor.py -x '+args.gfa + ' -l '+args.output+'/new_links_scaled_iteration_'+str(iter_num)+'_sorted -f '+ args.length+' -c '+args.cutoff+' -i '+str(iter_num)+ ' -s '+args.length + ' -g '+args.bed + ' -b '+args.output+'/iteration_'+str(iter_num+1)+'.bed -n '+args.output+'/scaffolds_'+str(iter_num)+'.p'
			log.write(cmd+'\n')
			p = subprocess.check_output(cmd,shell=True)

		except subprocess.CalledProcessError as err:
			print >> sys.stderr, str(err.output)
			sys.exit(1)

	#Now reconcile all the information to get new RE counts and lengths
	if not os.path.isfile(args.output+'/scaffold_length_iteration_'+str(iter_num+1)):
		try:
			cmd = 'python reconcile.py -p '+args.output+'/scaffolds_'+str(iter_num)+'.p -l '+args.length + ' -i '+str(iter_num+1) + ' -r '+args.output+'/RE_sites_iteration_1 -o '+args.output+'/scaffold_length_iteration_'+str(iter_num+1)+' -x '+args.output+'/RE_sites_iteration_'+str(iter_num+1)
			log.write(cmd+'\n')
			p = subprocess.check_output(cmd,shell=True)
		except subprocess.CalledProcessError as err:
			print >> sys.stderr, str(err.output)
			sys.exit(1)

	iter_num += 1

	#now do iterative
	while True:
		if not os.path.isfile(args.output+'/new_links_iteration_'+str(iter_num)):
			try:
				cmd = 'python new_links.py -b '+ args.output+'/iteration_'+str(iter_num)+'.bed' + ' -c '+ args.output+'/RE_sites_iteration_'+str(iter_num)+ ' -l '+args.output+'/scaffold_length_iteration_'+str(iter_num)+ ' > '+args.output+'/new_links_iteration_'+str(iter_num)
				log.write(cmd+'\n')
				p = subprocess.check_output(cmd,shell=True)

			except subprocess.CalledProcessError as err:
				print >> sys.stderr, str(err.output)
				sys.exit(1)


		if not os.path.isfile(args.output+'/new_links_scaled_iteration_'+str(iter_num)):
			try:
				cmd = 'python scaled_scores.py '+args.output+'/new_links_iteration_'+str(iter_num)+ ' > ' + args.output+'/new_links_scaled_iteration_'+str(iter_num)
				log.write(cmd+'\n')
				p = subprocess.check_output(cmd,shell=True)
			except subprocess.CalledProcessError as err:
				print >> sys.stderr, str(err.output)
				sys.exit(1)

		#Sort the links by column 5 
		if not os.path.isfile(args.output+'/new_links_scaled_iteration_'+str(iter_num)+'_sorted'):
			try:
				cmd = 'sort -k 5 -gr '+args.output+'/new_links_scaled_iteration_'+str(iter_num)+ ' > '+ args.output+'/new_links_scaled_iteration_'+str(iter_num)+'_sorted'
				log.write(cmd+'\n')
				p = subprocess.check_output(cmd,shell=True)
			except subprocess.CalledProcessError as err:
				print >> sys.stderr, str(err.output)
				sys.exit(1)


		#NOW check if any useful link here
		if check(args.output+'/new_links_scaled_iteration_'+str(iter_num)+'_sorted'):
			break

		if not os.path.isfile(args.output+'/scaffolds_'+str(iter_num)+'.p'):
			try:
				cmd = 'python layout_refactor.py -p '+args.output+'/scaffolds_'+str(iter_num-1)+'.p -f '+args.length + ' -l '+args.output+'/new_links_scaled_iteration_'+str(iter_num)+'_sorted -c '+args.cutoff+' -i '+str(iter_num)+ ' -s '+args.output+'/scaffold_length_iteration_'+str(iter_num) + ' -g '+args.bed + ' -b '+args.output+'/iteration_'+str(iter_num+1)+'.bed -n '+args.output+'/scaffolds_'+str(iter_num)+'.p'
				log.write(cmd+'\n')
				p = subprocess.check_output(cmd,shell=True)

			except subprocess.CalledProcessError as err:
				print >> sys.stderr, str(err.output)
				sys.exit(1)

		if not os.path.isfile(args.output+'/RE_sites_iteration_'+str(iter_num+1)):
			try:
				cmd = 'python reconcile.py -p '+args.output+'/scaffolds_'+str(iter_num)+'.p -l '+args.length + ' -i '+str(iter_num+1) + ' -r '+args.output+'/RE_sites_iteration_1 -o '+args.output+'/scaffold_length_iteration_'+str(iter_num+1)+' -x '+args.output+'/RE_sites_iteration_'+str(iter_num+1)
				log.write(cmd+'\n')
				p = subprocess.check_output(cmd,shell=True)
			except subprocess.CalledProcessError as err:
				print >> sys.stderr, str(err.output)
				sys.exit(1)

		iter_num += 1

	log.close()

if __name__ == '__main__':
	main()

