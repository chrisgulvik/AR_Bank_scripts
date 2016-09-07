#!/usr/bin/env python

__version__ = '1.0'


import argparse
import errno
import gzip
import logging
import os
import pwd
import re
import shutil
import subprocess
import sys
import textwrap as _textwrap
import threading
from datetime import datetime
from glob import glob
from tempfile import gettempdir
from multiprocessing import Process
from uuid import uuid4

def parseArgs():
	parser = argparse.ArgumentParser(add_help=False,
		formatter_class=ModifiedHelpFormatter,
		description='Performs magic, in parallel',
		epilog='''Output directory structure: L|
		1.QA_Reads L| 1.Reads L| 2.Kraken L|
		3.Asm L| 3.QA_Asm L|
		4.16S L| 4.AR L| 4.MLST L| 4.Plasmid L|
		5.Viz ''')
	req = parser.add_argument_group('Required')
	req.add_argument('-1', '--R1', required=True,
		help='input R1 fastq[.gz] file')
	req.add_argument('-2', '--R2', required=True,
		help='input R2 fastq[.gz] file')
	opt = parser.add_argument_group('Optional')
	opt.add_argument('-b', '--base', default=None,
		help='N|string to append to output files\n'
		'  default: first item split by underscore of R1 and R2 basenames')
	opt.add_argument('-c', '--cpus', type=str, default='2',
		help='N|number of CPUs\n'
		'  default: 2')
	opt.add_argument('-h', '--help', action='help',
		help='show this help message and exit')
	opt.add_argument('-o', '--outpath', default=None,
		help='N|output folder\n'
		'  default: `cwd`/basename/')
	opt.add_argument('-v', '--version', action='version',
		version='%(prog)s v{}'.format(__version__))
	opt.add_argument('--db_16S', default=None,
		help='N|FastA file of a 16S sequence\n'
		'  default: first file in exported $AR_BANK_DBs PATH(s) matching: *16S*.f*')
	opt.add_argument('--db_adapters', default=None,
		help='N|FastA file of adapters\n'
		'  default: first file in exported $AR_BANK_DBs PATH(s) matching: *dapter*.f*')
	opt.add_argument('--db_ar', default=None,
		help='N|FastA file of srst2-formatted AR genes\n'
		'  default: first file in exported $AR_BANK_DBs PATH(s) matching: *ARG-ANNOT*.f*')
	opt.add_argument('--db_kraken', default=None,
		help='N|kraken database (path) to use\n'
		'  default: first path in exported $KRAKEN_DEFAULT_DB or $KRAKEN_DB_PATH')
	opt.add_argument('--db_phiX', default=None,
		help='N|FastA file of PhiX genome\n'
		'  default: first file in exported $AR_BANK_DBs PATH(s) matching: *hiX*.f*')
	opt.add_argument('--db_plasmid', default=None,
		help='N|FastA file of PlasmidFinder sequences\n'
		'  default: first file in exported $AR_BANK_DBs PATH(s) matching: *lasmid*.f*')
	return parser.parse_args()

###############################################################################
###############################################################################

class ModifiedHelpFormatter(argparse.HelpFormatter):
	def _split_lines(self, text, width):
		''' interprets escape characters like newline and tab,
		when given within an add_argument help of ArgumentParser '''
		if text.startswith('N|'):
			return text[2:].splitlines()  
		return argparse.HelpFormatter._split_lines(self, text, width)

	def _fill_text(self, text, width, indent):
		''' creates newline+space*4 in description and epilog of ArgumentParser '''
		text = self._whitespace_matcher.sub(' ', text).strip()
		paragraphs = text.split('L|')
		multiline_text = ''
		for p in paragraphs:
			formatted_paragraph = _textwrap.fill(p, width,
				initial_indent=indent, subsequent_indent=indent) + '\n    '
			multiline_text = multiline_text + formatted_paragraph
		return multiline_text

class LogPipe(threading.Thread):
	def __init__(self, level):
		''' sets up the object with a logger and starts up the thread '''
		threading.Thread.__init__(self)
		self.daemon = False
		self.level = level
		self.fd_read, self.fd_write = os.pipe()
		self.pipe_reader = os.fdopen(self.fd_read)
		self.start()

	def fileno(self):
		''' returns the write file descriptor of the pipe '''
		return self.fd_write

	def run(self):
		''' runs the thread and logs everything '''
		for line in iter(self.pipe_reader.readline, ''):
			# logging.debug(line.strip('\n'))
			logging.log(self.level, line.strip('\n'))
		self.pipe_reader.close()

	def close(self):
		''' closes the write end of the pipe '''
		os.close(self.fd_write)

###############################################################################
###############################################################################

def dependency(dep):
	''' checks for binary or script availability and
	returns the path if it is available '''
	for path in os.environ.get('PATH', '').split(':'):
		if os.path.exists(os.path.join(path, dep)) and \
		not os.path.isdir(os.path.join(path, dep)):
			return os.path.join(path, dep)
	return None

def cat_it(infiles, outfile):
	''' combines the contents of a list of files into a single output file '''	
	with open(outfile, 'w') as o:
		for f in infiles:
			shutil.copyfileobj(open(f, 'r'), o)

def verify_file(testfile):
	''' verifies file is present and non-zero size '''
	if not os.path.isfile(testfile):
		logging.error('{} file absent'.format(testfile))
		sys.exit('ERROR: {} file absent'.format(testfile))
	if os.path.getsize(testfile) == 0:
		logging.error('{} file present but lacks content'.format(testfile))
		sys.exit('ERROR: {} file present but lacks content'.format(testfile))

def syscall(syscmd, std_in, out, err, cwd_path):
	''' given a syscmd list, executes a system call 
	and optionally logs stdout and stderr '''
	if out == 'dump' or err == 'dump':
		with open(os.devnull, 'wb') as dump:
			if out == 'dump' and err == 'dump':
				returncode = subprocess.call(syscmd,
					stdin=std_in, stdout=dump, stderr=dump, cwd=cwd_path, shell=False)
				if returncode != 0:
					logging.error('failed system call ' + ' '.join(syscmd))
					sys.exit('ERROR: failed system call\ncheck the logfile')
				logging.info('completed system call ' + ' '.join(syscmd))
			elif out == 'dump':
				e = LogPipe(logging.DEBUG)
				s = subprocess.Popen(syscmd,
					stdin=std_in, stdout=dump, stderr=e, cwd=cwd_path, shell=False)
				if s.wait() != 0:  #return code
					logging.error('failed system call ' + ' '.join(syscmd))
					e.close()
					sys.exit('ERROR: failed system call\ncheck the logfile')
				e.close()
				logging.info('completed system call ' + ' '.join(syscmd))
			elif err == 'dump':
				o = LogPipe(logging.DEBUG)
				s = subprocess.Popen(syscmd,
					stdin=std_in, stdout=o, stderr=dump, cwd=cwd_path, shell=False)
				if s.wait() != 0:  #return code
					logging.error('failed system call ' + ' '.join(syscmd))
					o.close()
					sys.exit('ERROR: failed system call\ncheck the logfile')
				o.close()
				logging.info('completed system call ' + ' '.join(syscmd))
	else:
		o = LogPipe(logging.DEBUG)
		e = LogPipe(logging.DEBUG)
		s = subprocess.Popen(syscmd,
			stdin=std_in, stdout=o, stderr=e, cwd=cwd_path, shell=False)
		if s.wait() != 0:  #return code
			logging.error('failed system call ' + ' '.join(syscmd))
			o.close()
			e.close()
			sys.exit('ERROR: failed system call\ncheck the logfile')
		o.close()
		e.close()
		logging.info('completed system call ' + ' '.join(syscmd))

def shell_syscall(shell_cmd, std_in, out, err, cwd_path):
	''' given a shell_cmd string, executes system calls
	where shell is required '''
	if std_in == 'pipe':
		std_in = subprocess.PIPE
	if out == 'pipe':
		out = subprocess.PIPE
	if err == 'pipe':
		err = subprocess.PIPE

	if out == 'dump' or err == 'dump':
		with open(os.devnull, 'wb') as dump:
			if out == 'dump' and err == 'dump':
				r = subprocess.Popen(shell_cmd,
					stdin=std_in, stdout=dump, stderr=dump, cwd=cwd_path, shell=True)
			elif out == 'dump':
				r = subprocess.Popen(shell_cmd,
					stdin=std_in, stdout=dump, stderr=err, cwd=cwd_path, shell=True)
			elif err == 'dump':
				r = subprocess.Popen(shell_cmd,
					stdin=std_in, stdout=out, stderr=dump, cwd=cwd_path, shell=True)
	else:
		r = subprocess.Popen(shell_cmd,
			stdin=std_in, stdout=out, stderr=err, cwd=cwd_path, shell=True)

	if r.wait() == 0:
		logging.info('completed system call ' + shell_cmd)
	else:
		print 'RETURN CODE={}'.format(r)
		logging.error('failed system call ' + shell_cmd)
		sys.exit('ERROR: failed system call\ncheck the logfile')

def bytes_to_human_readable(num_bytes):
	''' given a int of bytes, returns a human-readable string based on:
	1 kilobyte (kB) = 2^10 bytes = 1,024 B
	1 megabyte (MB) = 2^20 bytes = 1,048,576 B
	1 gigabyte (GB) = 2^30 bytes = 1,073,741,824 B '''
	if num_bytes == 0:
		return '0 B'
	units = ['B', 'KB', 'MB', 'GB']
	i = 0
	while num_bytes >= 1024 and i < len(units)-1:
		num_bytes /= 1024.
		i += 1
	val = ('{0:.1f}'.format(num_bytes)).rstrip('0').rstrip('.')
	return '{} {}'.format(val, units[i])

def gzip_it(infile, outfile):
	# ''' gunzip compresses a file (gzip -fk) '''
	# with open(infile, 'rb') as i_file:
	# 	with gzip.open(outfile, 'wb') as o_file:
	# 		shutil.copyfileobj(i_file, o_file)
	shell_syscall('pigz -f -9 {} > {}'.format(infile, outfile),
		None, 'pipe', 'dump', None)

###############################################################################
###############################################################################

def interpose_raw_reads(R1, R2, outfile):
	''' interleaves sister reads into a single FastQ output file '''
	# this function is the work of Nick Crawford and not mine          #
	# source1: http://gist.github.com/ngcrawford/2232505               #
	# source2: older version: 080e07f2b4759baf4831bb8e3cb67dc9ac158027 #
	# obtained: Fri 12 Feb 2016                                        #
	# merged sources 1 and 2; mostly stylistic modifications           #
	if R1.endswith('.gz') and R2.endswith('.gz'):
		left  = gzip.open(R1, 'rb')
		right = gzip.open(R2, 'rb')
	else:
		left  = open(R1, 'rU')
		right = open(R2, 'rU')
	fout  = open(outfile, 'wb')
	while 1:
		# process the R1 file
		left_id = left.readline()
		if not left_id: 
			break
		left_seq    = left.readline()
		left_plus   = left.readline()
		left_quals  = left.readline()
		# process the R2 file
		right_id    = right.readline()
		right_seq   = right.readline()
		right_plus  = right.readline()
		right_quals = right.readline()
		# write output
		fout.write(left_id)
		fout.write(left_seq)
		fout.write(left_plus)
		fout.write(left_quals)
		fout.write(right_id)
		fout.write(right_seq)
		fout.write(right_plus)
		fout.write(right_quals)
	left.close()
	right.close()
	fout.close()
	logging.info('Interposing reads for NCBI SRA archiving complete')

def run_QA_raw_reads(QA_raw_script, R1, R2, QA_seqs_outdir):
	''' calculates read quality metrics from raw Illumina 1.8 
	paired-end reads '''
	cmd = '{} {} {} {}'.format(QA_raw_script, R1, R2, QA_seqs_outdir)
	shell_syscall(cmd, None, 'pipe', 'dump', None)
	logging.info('QA of raw reads complete')

def run_bbduk(bbduk, R1, R2, phiX, trim_outpath, base, cpus):
	''' runs Brian Bushnell's Decontamination Using Kmers (BBDUK) to 
	remove PhiX from reads '''
	noPhiX_R1 = os.path.join(trim_outpath, base + '-noPhiX-R1.fsq')
	noPhiX_R2 = os.path.join(trim_outpath, base + '-noPhiX-R2.fsq')
	cmd = [bbduk, 'threads=' + cpus, 'in=' + R1, 'in2=' + R2,
			'out=' + noPhiX_R1, 'out2=' + noPhiX_R2,
			'ref=' + phiX, 'k=31', 'hdist=1', 'overwrite=t']
	syscall(cmd, None, 'log', 'log', None)
	logging.info('BBDuk complete')
	return noPhiX_R1, noPhiX_R2

def run_trimmo(trimmo, noPhiX_R1, noPhiX_R2, adapters, trim_outpath, base,
	cpus):
	''' clips off adapters, performs a sliding window quality trim,
	and discards short remaining reads '''
	cmd = ('{} PE -phred33 -threads {} {} {} {} {} {} {} '
			'ILLUMINACLIP:{}:2:20:10:8:TRUE SLIDINGWINDOW:20:30 LEADING:20 '
			'TRAILING:20 MINLEN:50').format(trimmo, cpus, noPhiX_R1, noPhiX_R2,
			os.path.join(trim_outpath, base + '_R1.paired.fq'),
			os.path.join(trim_outpath, base + '_R1.unpaired.fq'),
			os.path.join(trim_outpath, base + '_R2.paired.fq'),
			os.path.join(trim_outpath, base + '_R2.unpaired.fq'), adapters)
	syscall(cmd.split(), None, 'dump', 'log', None)
	unpaired_files = [os.path.join(trim_outpath, base + '_R1.unpaired.fq'),
		os.path.join(trim_outpath, base + '_R2.unpaired.fq')]
	cat_it(unpaired_files, os.path.join(trim_outpath, base + '_single.fq'))
	os.remove(os.path.join(trim_outpath, base + '_R1.unpaired.fq'))
	os.remove(os.path.join(trim_outpath, base + '_R2.unpaired.fq'))
	os.remove(os.path.join(trim_outpath, base + '-noPhiX-R1.fsq'))
	os.remove(os.path.join(trim_outpath, base + '-noPhiX-R2.fsq'))
	clean_R1         = os.path.join(trim_outpath, base + '_R1.paired.fq')
	clean_R2         = os.path.join(trim_outpath, base + '_R2.paired.fq')
	clean_singletons = os.path.join(trim_outpath, base + '_single.fq')
	logging.info('Trimmomatic complete')
	return clean_R1, clean_R2, clean_singletons

def run_QA_clean_reads(QA_clean_script, R1, R2, single, QA_seqs_outdir):
	''' calculates read quality metrics from cleaned Illumina 1.8
	paired-end reads with a singleton read file as well '''
	cmd = '{} {} {} {} {}'.format(QA_clean_script, R1, R2, single,
			QA_seqs_outdir)
	shell_syscall(cmd, None, 'pipe', 'dump', None)
	logging.info('QA of cleaned reads complete')

def run_kraken(kraken, kraken_report, sum_script, clean_R1, clean_R2, db,
	kraken_outpref, cpus, tmp_path):
	''' executes kraken and summarizes the report to assess contamination by
	listing read quantities that match to unclassified, top 3 species, 
	and top 3 genera '''
	kraken_outfile = kraken_outpref + '.full_report.tab'
	cmd1 = ('{} --threads {} --db {} --preload --paired --fastq-input '
			'{} {} | {} > {}').format(kraken, cpus, db, clean_R1,
			clean_R2, kraken_report, kraken_outfile)
	shell_syscall(cmd1, None, 'pipe', 'dump', None)
	summary_file = kraken_outpref + '.Summary.tab'
	cmd2 = '{} {} > {}'.format(sum_script, kraken_outfile, summary_file)
	shell_syscall(cmd2, None, 'pipe', 'dump', None)
	logging.info('Kraken complete')

def run_spades(spades, clean_R1, clean_R2, clean_singletons, spades_outpath,
	cpus, tmp_path):
	''' generates an assembly (without kmer optimization) '''
	cmd = [spades, '--pe1-1', clean_R1, '--pe1-2', clean_R2,
			'--pe1-s', clean_singletons, '-o', spades_outpath,
			'--phred-offset', '33', '--only-assembler', #'-k', '29,55,97',
			'--cov-cutoff', 'auto', '-t' + cpus]
	syscall(cmd, None, 'dump', 'log', None)
	Kmer_dirs = glob(os.path.join(spades_outpath, 'K*'))
	for d in Kmer_dirs:
		shutil.rmtree(d)
	shutil.rmtree(os.path.join(spades_outpath, 'misc'))
	shutil.rmtree(os.path.join(spades_outpath, 'tmp'))
	for junk in ['before_rr.fasta', 'contigs.paths', 'dataset.info',
		'input_dataset.yaml', 'params.txt', 'scaffolds.fasta',
		'scaffolds.paths']:
		os.remove(os.path.join(spades_outpath, junk))
	logging.info('SPAdes complete')

def run_QA_asm(quast, asm, quast_outpath, cpus):
	''' calculates several metrics from raw assembly,
	such as:  N50, # genes, # contigs, and cumulative length '''
	cmd = [quast, '-o', quast_outpath, asm, '--min-contig', '500',
		'--no-html', '--no-snps', '--gene-finding',
		'--ambiguity-usage', 'one', '--threads', cpus]
	syscall(cmd, None, 'log', 'log', None)
	for junk in ['report.tex', 'transposed_report.tex', 'report.txt',
		'transposed_report.txt']:
		os.remove(os.path.join(quast_outpath, junk))
	shutil.rmtree(os.path.join(quast_outpath, 'predicted_genes'))
	os.rmdir(os.path.join(quast_outpath, 'basic_stats'))
	logging.info('Quast complete')

def clean_asm(filt_script, init_asm, clean_asm_outfile, clean_asm_outpath,
	base):
	''' filters out low quality contigs '''
	cmd = '{} -i {} -o {} -p {} -b {}'.format(filt_script,
	init_asm, clean_asm_outfile, clean_asm_outpath, base)
	shell_syscall(cmd, None, 'dump', 'dump', None)
	logging.info('Contig filtering complete')

def annot_asm(prokka, clean_asm, base, clean_asm_outpath, cpus, tmp_path):
	''' produces an annotated genbank file from an assembly '''
	cleaned_base = re.sub(r'[^\w]', '', base).replace('_', '')
	cmd = [prokka, '--outdir', tmp_path, '--force', '--addgenes', 
			'--locustag', cleaned_base, '--prefix', base,
			'--cpus', cpus, '--mincontiglen', '500', '--evalue', '1e-06', 
			'--rnammer', clean_asm]
	syscall(cmd, None, 'log', 'dump', None)
	for e in ['faa','ffn','fna','fsa','gff','log','sqn','tbl','txt']:
		os.remove(os.path.join(tmp_path, base + '.' + e))
	try:
		os.remove(os.path.join(tmp_path, base + '.err'))
	except OSError:
		pass
	shutil.move(os.path.join(tmp_path, base + '.gbk'), 
		os.path.join(clean_asm_outpath, base + '.gbk'))
	logging.info('Genome annotation complete')

def extract_16S(rnammer, blastn, asm, rnammer_outpref, cpus, tmp_path):
	''' extracts all 16S rRNA gene sequences from an assembly
	and uses ports 4444 and 4544 to connect to NCBI '''
	rnammer_outfile  = rnammer_outpref + '16S.fa'
	blastn_outfile   = rnammer_outpref + '16S.nt.blastn.tab'
	species_outfile  = rnammer_outpref + 'species.tab'
	cmd1 = '{} -S bac -m ssu -f {} < {}'.format(rnammer, rnammer_outfile, asm)
	shell_syscall(cmd1, 'pipe', None, None, tmp_path)
	verify_file(rnammer_outfile)
	fmt = ('6 qseqid sseqid pident length mismatch gapopen qstart qend '
			'sstart send evalue bitscore qlen')
	cmd2 = [blastn, '-word_size', '10', '-task', 'blastn', '-remote',
			'-db', 'nt', '-max_hsps', '1', '-max_target_seqs', '1',
			'-query', rnammer_outfile, '-out', blastn_outfile,
			'-outfmt', fmt]
	syscall(cmd2, None, 'log', 'log', tmp_path)
	verify_file(blastn_outfile)
	with open(blastn_outfile, 'r') as fin:
		top_hit = fin.readline().rstrip()
	species_name = top_hit.split('\t')[-1]
	hit_length   = top_hit.split('\t')[3]
	with open(species_outfile, 'w') as fout:
		fout.write('{}\t{}'.format(hit_length, species_name))
	logging.info('RNAmmer complete')

def run_cSSTAR(cSSTAR, asm, AR_db, base, cSSTAR_outpath, tmp_path):
	''' identifies antibiotic resistance determinants from a genome assembly '''
	cmd = '{} -g {} -d {} -b {}'.format(cSSTAR, asm, AR_db, base)
	with open(os.path.join(cSSTAR_outpath, base + '.ARG-ANNOT.tab'), 'w') as o:
		syscall(cmd.split(), None, o, 'dump', tmp_path)
	shutil.move(os.path.join(tmp_path, base + '.blastn.tsv'),
		os.path.join(cSSTAR_outpath, base + '.blastn.tsv'))
	os.remove(os.path.join(tmp_path, 'c-SSTAR_' + base + '.log'))
	logging.info('c-SSTAR complete')

def run_TSeeman_MLST(mlst, asm, asm_MLST_outfile, cpus):
	''' identifies MLST from an assembly without a species name given '''
	cmd = '{} --threads={} {} > {}'.format(mlst, cpus, asm, asm_MLST_outfile)
	shell_syscall(cmd, None, 'pipe', 'dump', None)
	logging.info('MLST complete')

def run_plasmidfinder(blastn, mkblastdb, asm, plasfind_outpref,
	plasfind_outfile, plasmid_db, cpus, tmp_path):
	cmd_mkdb  = [mkblastdb, '-in', asm, '-out', plasfind_outpref,
				'-dbtype', 'nucl']
	cmd_blast = [blastn, '-word_size', '11', '-task', 'blastn', 
				'-evalue', '1e-10', '-num_threads', cpus,
				'-query', plasmid_db, '-db', plasfind_outpref,
				'-out', plasfind_outfile, '-outfmt', '6']
	syscall(cmd_mkdb, None, 'log', 'log', tmp_path)
	syscall(cmd_blast, None, 'log', 'log', tmp_path)
	os.remove(plasfind_outpref + '.nhr')
	os.remove(plasfind_outpref + '.nin')
	os.remove(plasfind_outpref + '.nsq')
	logging.info('PlasmidFinder complete')

def run_bandage(bandage, graph, sekstenS_file, outpref):
	cmd1 = [bandage, 'image', graph, outpref + '_16S-in-Blue.png', '--query',
			sekstenS_file, '--scope', 'aroundblast', '--distance', '2']
	cmd2 = [bandage, 'image', graph, outpref + '_full.svg', '--fontsize', '10',
			'--colour', 'random', '--mindepth', '4', '--edgewidth', '5']
	syscall(cmd1, None, 'dump', 'log', None)
	syscall(cmd2, None, 'dump', 'log', None)
	logging.info('Bandage complete')

###############################################################################
###############################################################################

init_time = datetime.now()
opts = parseArgs()
base = opts.base
cpus = opts.cpus
R1, R2 = opts.R1, opts.R2

# I/O handling
if base is None:
	if os.path.basename(R1).split('_')[0] == \
	os.path.basename(R2).split('_')[0]:
		base = os.path.basename(R1).split('_')[0]
	else:
		sys.exit('ERROR: different basenames for R1 and R2 '
		' files split by first underscore')
inpath = os.path.dirname(os.path.abspath(R1))
if opts.outpath:
	outpath = os.path.abspath(opts.outpath)
elif opts.outpath is None:
	outpath = os.path.join(inpath, base)
if not os.path.exists(outpath):
	os.mkdir(outpath)

# Setup logging
logging.basicConfig(filename=os.path.join(outpath , base + '.log'),
	format='%(asctime)s: %(levelname)s: %(message)s', level=logging.DEBUG,
	filemode='w', datefmt='%a %d-%b-%Y %H:%M:%S')
logging.info('#### System information ####')
logging.info('user: {}'.format(pwd.getpwuid(os.getuid()).pw_name))
logging.info('release: {}'.format(os.uname()[3]))
logging.info('shell env: {}'.format(pwd.getpwuid(os.getuid()).pw_shell))
logging.info('cwd: {}'.format(os.getcwd()))
logging.info('python version: {}'.format(sys.version))

# Check binary and script dependencies and log paths
logging.info('#### Paths for binaries and scripts used ####')
lib = {}
for b in ['QualAssessRawSeqs.bash', 'bbduk.sh',
'trimmomatic', 'QualAssessCleanSeqs.bash', 'spades.py', 'kraken',
'kraken-report', 'summarize_kraken-report.sh', 'quast', 'filter.contigs.py',
'prokka', 'rnammer', 'c-SSTAR', 'mlst', 'blastn', 'makeblastdb', 'Bandage']:
	avail = dependency(b)
	if avail:
		dependency_path = os.path.realpath(avail)  #follow symlinks (important for trimmo but will just do it for all)
		logging.info('found {}'.format(dependency_path))
		lib[b] = dependency_path
	else:
		print '\tERROR: {} not found'.format(b)
		sys.exit(1)

# Get database files and log paths
logging.info('#### Paths for reference sequences and databases used ####')
d0 = [opts.db_phiX,     'db_phiX',     '*hiX*.f*',       'phiX_genome']
d1 = [opts.db_adapters, 'db_adapters', '*dapter*.f*',    'adapters']
d2 = [opts.db_ar,       'db_ar',       '*ARG-ANNOT*.f*', 'ar_db']
d3 = [opts.db_plasmid,  'db_plasmid',  '*lasmid*.f*',    'plasmid_db']
d4 = [opts.db_16S,      'db_16S',      '*16S*.f*',       'sekstenS_file']
FastA_DBs = {}
for d in [d0, d1, d2, d3, d4]:
	if d[0]:
		FastA_DBs[d[3]] = str(d[0])
	else:
		if os.environ.get('AR_BANK_DBs'):
			db_paths = os.environ.get('AR_BANK_DBs').split(':')
			finds = []
			for p in db_paths:
				search_found = glob(os.path.join(p, d[2]))
				if search_found:
					finds.append(search_found[0])
			if finds:
				FastA_DBs[d[3]] = str(finds[0])
			else:
				logging.error(('ERROR: AR_BANK_DBs export path lacks '
								'{} and it was not provided as an '
								'argument').format(d[1]))
				sys.exit(('ERROR: AR_BANK_DBs export path lacks {} '
						'and it was not provided as an arg').format(d[1]))
		else:
			logging.error(('ERROR: AR_BANK_DBs export path lacks '
							'{} and it was not provided as an '
							'argument').format(d[1]))
			sys.exit(('ERROR: AR_BANK_DBs export path lacks {} '
					'and it was not provided as an arg').format(d[1]))
	size = bytes_to_human_readable(os.path.getsize(FastA_DBs[d[3]])) #also verifies file exists and is accessible
	logging.info('using database: {} [{}]'.format(FastA_DBs[d[3]], size))
# adapters = os.path.expanduser('~/barcodes/adapters_Nextera_NEB_TruSeq_NuGEN.fas')
# phiX_genome = os.path.expanduser('~/genome_references/PhiX/PhiX_NC_001422.1.fasta')
# ar_db = os.path.expanduser('~/AR/ARG-ANNOT.srst2.fasta')
# plasmid_db = os.path.expanduser('~/AR/plasmid_database.fsa')
# sekstenS_file = os.path.expanduser('~/reference_genes/Ecoli_16S.fa')
# export AR_BANK_DBs="$HOME/genome_references/PhiX:$HOME/AR:$HOME/barcodes:$HOME/reference_genes"

# Log dependency versions
#to-do: trimmo and rnammer
logging.info('#### Versions of binaries and scripts used ####')
logging.info('QualAssessRawSeqs.bash ' + subprocess.check_output(
	'{} --version'.format(
	lib['QualAssessRawSeqs.bash']), shell=True).rstrip())
logging.info('QualAssessCleanSeqs.bash ' + subprocess.check_output(
	'{} --version'.format(
	lib['QualAssessCleanSeqs.bash']), shell=True).rstrip())
logging.info(subprocess.check_output(
	'{} --version | grep "version"'.format(
	lib['kraken']), shell=True).rstrip())
logging.info('Bandage ' + subprocess.check_output(
	'{} --version'.format(
	lib['Bandage']), shell=True).rstrip())
logging.info(subprocess.check_output(
	'{} --version 2>&1'.format(
	lib['c-SSTAR']), shell=True).rstrip())
logging.info(subprocess.check_output(
	'{} --version 2>&1'.format(
	lib['quast']), shell=True).rstrip())
logging.info(subprocess.check_output(
	'{} --version 2>&1'.format(
	lib['spades.py']), shell=True).rstrip())
logging.info(subprocess.check_output(
	'{} --version 2>&1'.format(
	lib['prokka']), shell=True).rstrip())
logging.info('BBDuk ' + subprocess.check_output(
	'{} --help | grep "Last modified"'.format(
	lib['bbduk.sh']), shell=True).rstrip())
logging.info(subprocess.check_output('{} --version'.format(
	lib['mlst']), shell=True).rstrip())
logging.info(subprocess.check_output('{} -version | grep "blastn:"\n'.format(
	lib['blastn']), shell=True).rstrip())

# Output directory structure
tmp_path = os.path.join(gettempdir(), str(uuid4()))
for d in ['1.QA_Reads', '1.Reads', '2.Kraken', '3.Asm', '3.QA_Asm',
'4.MLST', '4.AR', '4.16S', '4.Plasmid', '5.Viz', tmp_path]:
	try:
		os.mkdir(os.path.join(outpath, d))
	except OSError as e:
		if e.errno == errno.EEXIST:
			logging.warning('{} already exists'.format(
				os.path.join(outpath, d)))
			logging.warning('conflicting files will be overwritten')
		else:
			raise

# Use raw reads
QA_seqs_outdir = os.path.join(outpath, '1.QA_Reads')
raw_interposed = os.path.join(outpath, '1.Reads', base + '.fastq')
p0 = Process(target=run_QA_raw_reads,
	args=(lib['QualAssessRawSeqs.bash'], R1, R2, QA_seqs_outdir))
p1 = Process(target=interpose_raw_reads,
	args=(R1, R2, raw_interposed))
p0.start(), p1.start()
p0.join(), p1.join()
verify_file(raw_interposed)
gzip_it(raw_interposed, raw_interposed + '.gz')
verify_file(raw_interposed + '.gz')

# Remove PhiX, adapter clip, and quality trim; QA before and after
trim_outpath = os.path.join(outpath, '1.Reads')
noPhiX_R1, noPhiX_R2 = run_bbduk(lib['bbduk.sh'], R1, R2,
	FastA_DBs['phiX_genome'], trim_outpath, base, cpus)
clean_R1, clean_R2, clean_singletons = run_trimmo(lib['trimmomatic'],
	noPhiX_R1, noPhiX_R2, FastA_DBs['adapters'], trim_outpath, base, cpus)
run_QA_clean_reads(lib['QualAssessCleanSeqs.bash'], clean_R1, clean_R2,
	clean_singletons, QA_seqs_outdir)

# Run kraken
kraken_default_db = os.environ.get('KRAKEN_DEFAULT_DB')
kraken_db_path    = os.environ.get('KRAKEN_DB_PATH')
if opts.db_kraken:
	kraken_db = opts.db_kraken
elif kraken_default_db and os.path.isdir(kraken_default_db):
	kraken_db = kraken_default_db
elif kraken_db_path and os.path.isdir(kraken_db_path):
	kraken_db = kraken_db_path
else:
	logging.error('kraken database not exported in env nor given as an arg')
	sys.exit('ERROR: kraken database not exported in env nor given as an arg')
verify_file(os.path.join(kraken_db, 'database.kdb'))
kraken_outpref = os.path.join(outpath, '2.Kraken', base)
run_kraken(lib['kraken'], lib['kraken-report'], 
	lib['summarize_kraken-report.sh'], clean_R1, clean_R2,
	kraken_db, kraken_outpref, cpus, tmp_path)

# Assemble genome; QA before and after filtering
spades_outpath = os.path.join(outpath, '3.Asm', base + '.raw')
run_spades(lib['spades.py'], clean_R1, clean_R2, clean_singletons,
	spades_outpath, cpus, tmp_path)
init_asm = os.path.join(outpath, '3.Asm', base + '.raw', 'contigs.fasta')
verify_file(init_asm)
init_quast_outpath = os.path.join(outpath, '3.QA_Asm', base + '.init')
run_QA_asm(lib['quast'], init_asm, init_quast_outpath, cpus)
clean_asm_outfile = base + '.fna'
clean_asm_outpath = os.path.join(outpath, '3.Asm')
clean_asm(lib['filter.contigs.py'], init_asm, clean_asm_outfile,
	clean_asm_outpath, base)
cleaned_asm = os.path.join(outpath, '3.Asm', base + '.fna')
verify_file(cleaned_asm)
annot_asm(lib['prokka'], cleaned_asm, base, clean_asm_outpath, cpus,
	tmp_path)
clean_quast_outpath = os.path.join(outpath, '3.QA_Asm', base + '.clean')
run_QA_asm(lib['quast'], cleaned_asm, clean_quast_outpath, cpus)

# Extract info from genome (16S, MLST, plasmids)
rnammer_outpref = os.path.join(outpath, '4.16S', base + '.')
cSSTAR_outpath  = os.path.join(outpath, '4.AR')
p0 = Process(target=extract_16S,
	args=(lib['rnammer'], lib['blastn'], cleaned_asm, rnammer_outpref,
	cpus, tmp_path))
p1 = Process(target=run_cSSTAR,
	args=(lib['c-SSTAR'], cleaned_asm, FastA_DBs['ar_db'], base,
	cSSTAR_outpath, tmp_path))
p0.start(), p1.start()
p0.join(), p1.join()
verify_file(rnammer_outpref + '16S.fa')
asm_MLST_outfile = os.path.join(outpath, '4.MLST', base + '.mlst.tab')
run_TSeeman_MLST(lib['mlst'], cleaned_asm, asm_MLST_outfile, cpus)
plasfind_outpref = os.path.join(outpath, '4.Plasmid', base)
plasfind_outfile = plasfind_outpref + '.plasmid.blastn.tab'
run_plasmidfinder(lib['blastn'], lib['makeblastdb'], cleaned_asm,
	plasfind_outpref, plasfind_outfile, FastA_DBs['plasmid_db'], cpus,
	tmp_path)

# Create visualizations (show where 16S is in the de Bruijn graph)
graph   = os.path.join(outpath, '3.Asm', base + '.raw',
	'assembly_graph.fastg')
outpref = os.path.join(outpath, '5.Viz', base + '.asm_graph_')
run_bandage(lib['Bandage'], graph, FastA_DBs['sekstenS_file'], outpref)

# Summarize data in all analyses
shutil.rmtree(tmp_path, ignore_errors=True)
logging.info('total run time was: {!s}'.format((datetime.now()-init_time)))
