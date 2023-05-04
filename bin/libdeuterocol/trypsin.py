import os
import pathlib

class BrokenInterval(object):
	def __init__(self, spans=None):
		self.spans = [None, None] if spans is None else spans

	def __contains__(self, other):
		for span in self.spans:
			if span[0] is None and span[-1] is None: return True
			elif span[0] is None and other <= span[-1]: return True
			elif span[-1] is None and other >= span[0]: return True
			elif (span[0] is not None) and (span[-1] is not None) and (span[0] <= other <= span[-1]): return True
		return False

class Trypsin(object):
	def __init__(self, deut2):
		self.deut2 = deut2

	def cut_pdbs(self, outdir):
		if not os.path.isdir(outdir): os.mkdir(outdir)

		for fragment in self.deut2.get_fragments(): 
			#infn = '{tmdatadir}/pdb/{pdb}.pdb'.format(tmdatadir=self.deut2.tmdatadir, pdb=fragment['pdbc'][:4].lower())
			infn = self.deut2.pdbdir.joinpath(fragment['pdbc'][:4].lower() + '.pdb')
			outfn = '{outdir}/{pdbc}_{label}.pdb'.format(outdir=outdir, **fragment)
			with open(infn) as infh:
				with open(outfn, 'w') as outfh:
					cut_pdb(infh, fragment['span'], chain=fragment['pdbc'][-1], outfh=outfh)

def cut_pdb(pdbfh, spans, chain, outfh=None):
	allowed = BrokenInterval(spans)

	out = ''

	for l in pdbfh:
		if l.startswith('ATOM  ') or l.startswith('HETATM'):
			if chain == 'A' and l[21] == ' ': pass
			elif l[21] != chain: continue

			resi = int(l[22:26])
			if resi in allowed: out += l

	out += 'END\n'

	if outfh is not None: outfh.write(out)

	return out
