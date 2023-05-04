import urllib.error
import urllib.request
import os
import csv
import re
import subprocess
from Bio import PDB, SeqIO, Seq
import warnings
from kevtools_common.types import TCID
import json

VERBOSITY = 0

class OPMTMData(object):

	dbpath = None
	force = False

	def __init__(self, dbpath):
		self.dbpath = dbpath

		if not VERBOSITY:
			warnings.filterwarnings('ignore', '.')

	def fetch(self, force=None):
		if force is not None: self.force = force

		if not os.path.isdir(self.dbpath): os.mkdir(self.dbpath)

		for table in ('primary_structures', 'secondary_representations', 'structure_subunits'):
			outfn = '{}/{}.csv'.format(self.dbpath, table)

			if self.force or not os.path.isfile(outfn):
				with urllib.request.urlopen('https://lomize-group-opm.herokuapp.com//{}?fileFormat=csv'.format(table)) as infh:
					with open(outfn, 'wb') as outfh:
						outfh.write(infh.read())

		if self.force or not (os.path.isfile('{}/all_pdbs.tar.gz'.format(self.dbpath)) or (os.path.isdir('{}/pdb'.format(self.dbpath)))):
			with urllib.request.urlopen('https://storage.googleapis.com/opm-assets/pdb/tar_files/all_pdbs.tar.gz') as infh:
				with open('{}/all_pdbs.tar.gz'.format(self.dbpath), 'wb') as outfh:
					outfh.write(infh.read())

		if self.force or not os.path.isfile('secondaries.tsv'):
			self.collect_secondaries()

		if self.force or not os.path.isdir('{}/pdb'.format(self.dbpath)):
			os.mkdir('{}/pdb'.format(self.dbpath))
			subprocess.call(['tar', 'xzf', '{}/all_pdbs.tar.gz'.format(self.dbpath), '-C', self.dbpath])
		self.download_missing()

		if self.force or not os.path.isfile('{}/ASSIGNMENTS.TSV'.format(self.dbpath)):
			self.write_assignments()

		if self.force or not os.path.isfile('{}/sequences.faa'.format(self.dbpath)):
			self.write_sequences()

		if self.force or not os.path.isfile('{}/tcmap.tsv'.format(self.dbpath)):
			self.tcblast()

	def collect_secondaries(self):
		if not os.path.isdir('{}/detailed'.format(self.dbpath)): os.mkdir('{}/detailed'.format(self.dbpath))

		pdbpattern = re.compile('[0-9]...(?:_[A-Za-z])?')
		fetchme = {}
		with open('{}/primary_structures.csv'.format(self.dbpath)) as fh:
			keys = None
			for i, line in enumerate(csv.reader(fh)):
				if not i: keys = line

				else: 
					obj = dict(zip(keys, line))
					pdbid = re.findall(pdbpattern, obj['pdbid'])[0]
					opmid = obj['id']

					if int(obj['secondary_representations_count']) > 1: fetchme[pdbid] = opmid
		for pdbid in fetchme:
			if self.force or not os.path.isfile('{}/detailed/{}.json'.format(self.dbpath, pdbid)):
				with urllib.request.urlopen('http://lomize-group-opm.herokuapp.com//primary_structures/{}'.format(fetchme[pdbid])) as infh:
					with open('{}/detailed/{}.json'.format(self.dbpath, pdbid), 'w') as outfh:
						outfh.write(infh.read().decode('utf-8'))

		if self.force or not os.path.isfile('{}/secondaries.tsv'.format(self.dbpath)):
			secondaries = {}
			for pdbid in sorted(fetchme):
				with open('{}/detailed/{}.json'.format(self.dbpath, pdbid)) as infh:
					obj = json.load(infh)
					secondaries[pdbid.upper()] = [link['pdbid'].upper() for link in obj['secondary_representations']]
			with open('{}/secondaries.tsv'.format(self.dbpath), 'w') as fh:
				for pdbid in secondaries:
					fh.write('{}\t{}\n'.format(pdbid, ','.join(secondaries[pdbid])))

	def download_missing(self):
		with open('{}/secondaries.tsv'.format(self.dbpath)) as fh:
			secondaries = parse_ctsv(fh)

		allpdbs = set()
		for pdbid in secondaries:
			allpdbs.add(pdbid)
			for secpdbid in secondaries[pdbid]: allpdbs.add(secpdbid)

		pdbpattern = re.compile('[0-9]...(?:_[A-Za-z])?')
		with open('{}/primary_structures.csv'.format(self.dbpath)) as fh:
			keys = None
			for i, line in enumerate(csv.reader(fh)):
				if not i: keys = line
				else:
					obj = dict(zip(keys, line))
					pdbid = re.findall(pdbpattern, obj['pdbid'])[0].upper()
					allpdbs.add(pdbid)

		allpdbs = sorted(allpdbs)
		

		fetchme = []
		for pdbid in allpdbs:
			if not os.path.isfile('{}/pdb/{}.pdb'.format(self.dbpath, pdbid.lower())):
				fetchme.append(pdbid)


		removeme = []

		for pdbid in fetchme:
			#print('http://files.rcsb.org/download/{}.pdb.gz'.format(pdbid.lower()))
			if self.force or not os.path.isfile('{}/pdb/{}.pdb.gz'.format(self.dbpath, pdbid.lower())):
				try:
					with urllib.request.urlopen('http://files.rcsb.org/download/{}.pdb.gz'.format(pdbid.lower())) as infh:
						with open('{}/pdb/{}.pdb.gz'.format(self.dbpath, pdbid.lower()), 'wb') as outfh: 
							outfh.write(infh.read())
				except urllib.error.HTTPError:
					try: 
						with urllib.request.urlopen('http://files.rcsb.org/download/{}.cif.gz'.format(pdbid.lower())) as infh:
							with open('{}/pdb/{}.cif.gz'.format(self.dbpath, pdbid.lower()), 'wb') as outfh: 
								outfh.write(infh.read())
					except urllib.error.HTTPError: removeme.append(pdbid)
					
		for pdbid in removeme: fetchme.remove(pdbid)

		for pdbid in fetchme:
			if os.path.isfile('{}/pdb/{}.cif.gz'.format(self.dbpath, pdbid.lower())):
				subprocess.call(['gunzip', '{}/pdb/{}.cif.gz'.format(self.dbpath, pdbid.lower())])
				p = subprocess.Popen(['pdbset', 
					'xyzin', '{}/pdb/{}.cif'.format(self.dbpath, pdbid.lower()), 
					'xyzout', '{}/pdb/{}.pdb'.format(self.dbpath, pdbid.lower())], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
				p.communicate(input='WRITE PDB'.encode('utf-8'))
				
			else:
				subprocess.call(['gunzip', '{}/pdb/{}.pdb.gz'.format(self.dbpath, pdbid.lower())])

	def write_assignments(self):
		pdbpattern = re.compile('[0-9]...(?:_[A-Za-z])?')
		assignments = {}

		with open('{}/secondaries.tsv'.format(self.dbpath)) as fh:
			secondaries = parse_ctsv(fh)

		with open('{}/structure_subunits.csv'.format(self.dbpath)) as fh:
			keys = None
			for i, line in enumerate(csv.reader(fh)):
				if not i: 
					keys = line
				else:
					obj = dict(zip(keys, line))
					pdbid = re.findall(pdbpattern, obj['pdbid'])[0].upper()
					chain = obj['protein_letter']
					pdbc = '{}_{}'.format(pdbid.upper(), chain)

					segment = re.sub('[^-() 0-9,A-Z]', '', obj['segment'])

					rawindices = [x for x in re.findall('[0-9]+', segment)]
					if (len(rawindices) % 3): raise ValueError('Something is very wrong with .segment for {}'.format(pdbc))

					tmss = {}
					current = None
					for i, x in enumerate(rawindices):
						if (i % 3) == 0:
							current = int(x)
						elif (i % 3) == 1:
							if current not in tmss: tmss[current] = [int(x), None]
						else:
							tmss[current][1] = int(x)
					indices = [tmss[i] for i in sorted(tmss)]
					assignments[pdbc] = indices

					if pdbid in secondaries:
						for secpdb in secondaries[pdbid]:
							secpdbc = secpdb + '_{}'.format(chain)
							assignments[secpdbc] = indices

		with open('{}/ASSIGNMENTS.TSV'.format(self.dbpath), 'w') as fh:
			for pdbc in sorted(assignments):
				fh.write('{}\t{}\n'.format(pdbc, ','.join(['{}-{}'.format(*span) for span in assignments[pdbc]])))

	def write_sequences(self):
		sequences = {}
		for bn in sorted(os.listdir('{}/pdb'.format(self.dbpath))):
			if not bn.endswith('pdb'): continue
			fn = '{}/pdb/{}'.format(self.dbpath, bn)

			pdbid = bn[:4].upper()
			with open(fn) as fh:
				seq = parse_pdb_atom(fh)
			for chain in seq:
				sequences['{}_{}'.format(pdbid, chain)] = seq[chain]
			

		records = {}
		for pdbc in sequences:
			if len(sequences[pdbc]) < 10: continue
			if re.match('^X+$', sequences[pdbc]): continue
			records[pdbc] = SeqIO.SeqRecord(Seq.Seq(sequences[pdbc]), id=pdbc, description='')
		writeme = [records[pdbc] for pdbc in sorted(records)]
		SeqIO.write(writeme, '{}/sequences.faa'.format(self.dbpath), 'fasta')

	def write_sequences_pdbparser(self):
		sequences = {}
		parser = PDB.PDBParser()
		for bn in sorted(os.listdir('{}/pdb'.format(self.dbpath))):
			print(bn)
			fn = '{}/pdb/{}'.format(self.dbpath, bn)
			try: struc = parser.get_structure('currentstruc', fn)
			except Exception as e:
				print(fn)
				continue
			pdbid = bn[:4]
			seq = {}
			for model in struc:
				for chain in model:
					pdbc = '{}_{}'.format(pdbid, chain.id)
					seq[pdbc] = ''
					for residue in chain:
						try: seq[pdbc] += PDB.protein_letters_3to1[residue.get_resname()]
						except KeyError: continue
			for pdbc in seq:
				if len(seq[pdbc]) < 10: continue
				if re.match('^X+$', seq[pdbc]): continue
				sequences[pdbc] = SeqIO.SeqRecord(Seq.Seq(seq[pdbc]), id=pdbc, description='')
			del struc
		writeme = [sequences[pdbc] for pdbc in sorted(sequences)]
		SeqIO.write(writeme, '{}/sequences.faa'.format(self.dbpath), 'fasta')
			
	def tcblast(self):
		if self.force or not os.path.isfile('{}/tcdb.psq'.format(self.dbpath)):
			subprocess.call(['extractFamily.pl', '-i', 'all', '-f', 'blast', '-o', self.dbpath])

		if self.force or not os.path.isfile('{}/tcblast.tsv'.format(self.dbpath)): #TODO: Check for completeness too
			subprocess.call(['blastp', '-db', '{}/tcdb'.format(self.dbpath), '-query', '{}/sequences.faa'.format(self.dbpath), '-evalue', '1e-10', '-max_target_seqs', '10', '-outfmt', '6 std qseq sseq', '-out', '{}/tcblast.tsv'.format(self.dbpath), '-comp_based_stats', 'no', '-seg', 'no', '-num_threads', '4'])

		if self.force or not (os.path.isfile('{}/tcmap.tsv'.format(self.dbpath)) and os.path.isfile('{}/tcmap_beyond.tsv'.format(self.dbpath))):
			tcmap = {}
			tcmap_beyond = {}
			done = set()
			donepdb = set()
			with open('{}/tcblast.tsv'.format(self.dbpath)) as fh:
				for l in fh:
					pdbc, tcid, pident = l.split()[:3]
					expect = float(l.split()[10])
					qseq, sseq = l.split()[-2:]

					tcid = TCID(tcid)
					pident = float(pident)
					#if pdbc in tcmap.get(tcid, []): continue
					#elif pdbc in tcmap_beyond.get(tcid, []): continue
					if pdbc in done: continue
					elif expect <= 1e-50:
						try: tcmap[tcid].append(pdbc)
						except KeyError: tcmap[tcid] = [pdbc]
					elif (pident/100 * len(sseq.replace('-', ''))) < 20:
						try: tcmap_beyond[tcid].append(pdbc)
						except KeyError: tcmap_beyond[tcid] = [pdbc]
					elif pident < 65:
						try: tcmap_beyond[tcid].append(pdbc)
						except KeyError: tcmap_beyond[tcid] = [pdbc]
					else:
						try: tcmap[tcid].append(pdbc)
						except KeyError: tcmap[tcid] = [pdbc]
					done.add(pdbc)
					donepdb.add(pdbc[:4])

			with open('{}/secondaries.tsv'.format(self.dbpath)) as fh:
				secondaries = parse_ctsv(fh)
				for primary in secondaries:
					if primary not in donepdb:
						try: tcmap_beyond['0.X.0.0.0-X00000'].append(primary)
						except KeyError: tcmap_beyond['0.X.0.0.0-X00000'] = [primary]
						tcmap_beyond['0.X.0.0.0-X00000'].extend(secondaries[primary])

			pdbpattern = re.compile('[0-9]...(?:_[A-Za-z])?')
			with open('{}/primary_structures.csv'.format(self.dbpath)) as fh:
				keys = None
				for i, line in enumerate(csv.reader(fh)):
					if not i: keys = line
					else:
						obj = dict(zip(keys, line))
						pdbid = re.findall(pdbpattern, obj['pdbid'])[0].upper()
						if pdbid not in donepdb:
							if '0.X.0.0.0-X00000' not in tcmap_beyond:
								tcmap_beyond['0.X.0.0.0-X00000'] = [pdbid]
							elif pdbid not in tcmap_beyond['0.X.0.0.0-X00000']: 
								tcmap_beyond['0.X.0.0.0-X00000'].append(pdbid)

			with open('{}/tcmap.tsv'.format(self.dbpath), 'w') as fh:
				for i, tcid in enumerate(sorted(tcmap)):
					fh.write('{}\t{}\n'.format(tcid, ','.join(sorted(tcmap[tcid]))))
					if not (i % 100): fh.flush()
			with open('{}/tcmap_beyond.tsv'.format(self.dbpath), 'w') as fh:
				for i, tcid in enumerate(sorted(tcmap_beyond)):
					fh.write('{}\t{}\n'.format(tcid, ','.join(sorted(tcmap_beyond[tcid]))))
					if not (i % 100): fh.flush()

def parse_ctsv(fh):
	obj = {}
	for l in fh:
		k, v = l.strip().split('\t')
		obj[k] = v.split(',')
	return obj

def parse_pdb_atom(fh):
	seq = {}
	for l in fh:
		if not (l.startswith('ATOM') or l.startswith('HETATM')): continue
		aname = l[13:15]
		if aname != 'CA': continue
		resname = l[17:20]
		if resname not in PDB.protein_letters_3to1: continue
		chain = 'A' if l[21] == ' ' else l[21]

		if chain not in seq: seq[chain] = ''
		seq[chain] += PDB.protein_letters_3to1[resname]
	return seq
