from .types import TCID
import urllib.request
import re

class Superfamily(object):
	def __init__(self, tcset=None, name=None):
		self.tcset = set() if tcset is None else tcset
		self.name = "Untitled superfamily" if name is None else name

	def __len__(self):
		return len(self.tcset)

	def __contains__(self, value):
		if value in self.tcset: return True
		return False

	def add(self, value):
		""" add a TCID to this Superfamily """
		#if it's already added, skip
		if value in self.tcset: return

		#if it's lower than an existing TCID, skip
		for tcid in self.tcset:
			if value in tcid: return

		#if it's higher than an existing TCID, remove the old one
		for tcid in self.tcset:
			if tcid in value:
				self.tcset.remove(tcid)

		#add it
		self.tcset.add(value)

	def remove(self, value):
		""" remove a TCID from this Superfamily """
		#if it's in the tcset, just remove it and hope self.tcset hasn't been tampered with
		if value in self.tcset:
			self.tcset.remove(value)
			return

		#if it's higher than an existing TCID, remove it as well
		for tcid in self.tcset:
			if tcid in value: self.tcset.remove(tcid)

		#if it's lower than an existing TCID, there's not much that can be done right now
		#TODO: subtractive self.exclude set maybe for superfamilies with exceptions

	def __repr__(self):
		return "Superfamily(name={}, length={})".format(repr(self.name), len(self))

	def __str__(self):
		return repr(self)

class SuperfamilyCollection(object):
	def __init__(self):
		self.superfamilies = []

	def __contains__(self, value):
		""" test if a TCID is in any of the superfamilies in this collection """
		for superfamily in self.superfamilies:
			if value in superfamily: return True
		return False

	def __len__(self): return len(self.superfamilies)

	def __iter__(self): return iter(self.superfamilies)

	def add(self, value):
		""" add a Superfamily to this SuperfamilyCollection. Currently performs no validation """
		self.superfamilies.append(value)

	def remove(self, value):
		""" remove a Superfamily from this SuperfamilyCollection """
		self.superfamilies.remove(value)

	def get_superfamily(self, value):
		""" return the superfamily a TCID belongs to or a NoneType """
		for superfamily in self.superfamilies:
			if value in superfamily: return superfamily
		return None


def abbreviate(sfname):
	""" attempt to automatically abbreviate a superfamily name """
	sfname = sfname.replace(" Superfamily", "").replace(" superfamily", "")
	if "(" in sfname: return re.findall("\((.*)\)", sfname)[0]
	elif sfname.startswith("ABC"): return "ABC"
	return sfname

def fetch_superfamilies(level=3, shortname=False):
	""" fetch the most recent version of the TCDB superfamilies table as a SuperfamilyCollection """

	with urllib.request.urlopen("https://tcdb.org/cgi-bin/substrates/listSuperfamilies.py") as fh:
		return parse_superfamilies(fh, level=level, shortname=shortname)

def parse_superfamilies(s, level=3, shortname=False):
	""" parse a string containing the TCDB superfamilies table as a SuperfamilyCollection """

	out = SuperfamilyCollection()
	sfdict = {}
	for l in s.split("\n"):
		if l.startswith("#"): continue
		elif not l.strip(): continue
		sl = l.split("\t")
		tcid = TCID(sl[0])[:level]
		sfname = sl[4]
		if shortname:
			sfname = abbreviate(sfname)

		if sfname not in sfdict: 
			sfdict[sfname] = Superfamily(name=sfname)
			out.add(sfdict[sfname])

		sfdict[sfname].add(tcid)

	return out
