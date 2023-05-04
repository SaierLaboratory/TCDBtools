import subprocess
import sys
import numpy as np

def warn(*things):
	print('[WARNING]', *things, file=sys.stderr)

def error(*things):
	print('[ERROR]', *things, file=sys.stderr)
	exit(1)

class IndexedSequence(object):
	record = None
	strict = False

	def __init__(self, record=None, start=0):
		#where record is a SeqIO.SeqRecord object or something close enough
		self.record = record
		self.start = start

	def __len__(self): return len(self.record.seq)

	def get_self_map_to_gapless(self):
		return self.get_self_map(blacklist=set(('-', '.', '_', ' ', 'X')))

	def get_self_revmap_to_gapless(self):
		'''

		'''
		return self.get_self_revmap(blacklist=set(('-', '.', '_', ' ', 'X')))

	def get_self_map(self, whitelist=set(), blacklist=set()):
		''' Provide a mapping between string-level indices and sequence-level indices.

		Parameters
		----------
		whitelist : iterable, optional
			An iterable of characters to explicitly treat as residues. 
		blacklist : iterable, optional
			An iterable of characters to explicitly treat as gaps.

		See also
		--------
		get_self_revmap : Provides the reverse mapping
		'''
		index = self.start

		mapping = {}

		for pos, resn in enumerate(self.record.seq):
			mapping[pos] = index

			if resn in whitelist: index += 1
			elif resn in blacklist: pass
			else: index += 1

		return mapping

	def get_self_revmap(self, whitelist=set(), blacklist=set()):
		''' Provide a mapping between sequence-level indices and string-level indices

		Parameters
		----------
		whitelist : iterable, optional
			An iterable of characters to explicitly treat as residues. 
		blacklist : iterable, optional
			An iterable of characters to explicitly treat as gaps.

		See also
		--------
		get_self_revmap : Provides the forward mapping
		'''
		index = self.start
		mapping = {}

		for pos, resn in enumerate(self.record.seq):
			if resn in whitelist:
				mapping[index] = pos
				index += 1
			elif resn in blacklist: pass
			else:
				mapping[index] = pos
				index += 1
		return mapping

	def get_transitive_map(self, other, whitelist=set(), blacklist=set(), strict=None):
		strict = self.strict if strict is None else strict

		if len(self) != len(other): 
			if strict: raise IndexError('Sequence lengths do not match (ensure sequences are aligned or disable strict mode)')

		selfmap = self.get_self_map(whitelist, blacklist)
		othermap = other.get_self_map(whitelist, blacklist)

		mapping = {}
		for iself, iother in zip(sorted(selfmap), sorted(othermap)):
			mapping[iself] = othermap[iother]

		return mapping

def collapse(arr):
	''' Collapse a nan-laden array for operations that don't tolerate such things '''
	newarr = []
	gaps = []	
	for i, x in enumerate(arr):
		if np.isnan(x): gaps.append(i)
		else: newarr.append(x)

	return np.array(newarr), np.array(gaps)

def uncollapse(arr, gaps):
	''' Uncollapse a nan-stripped array after operations that don't tolerate such things '''
	if not len(gaps): return np.array(arr[:])
	else:
		newarr = []
		oldarr = list(arr)
		for i in range(len(arr) + len(gaps)):
			if i in gaps: newarr.append(np.nan)
			#FIXME: weird issue where oldarr runs empty somehow
			elif len(oldarr): newarr.append(oldarr.pop(0))
		return np.array(newarr)

def unroll_ranges(spanlist):
    out = []
    for start, end in spanlist: 
        out.extend(range(start, end+1))

    return out

def roll_ranges(indexlist):
    out = []
    for index in indexlist:
        if not out: out.append([index, index])
        elif index == (out[-1][-1] + 1): out[-1][-1] += 1
        else: out.append([index, index])
    return out
