from . import tmalign
from . import superpose
#from . import dali
#from . import cealign


def get_aligner(name):
	if name == 'tmalign': return tmalign.TMAlign
	elif name == 'superpose' or name == 'ssm': return superpose.Superpose

	else: raise ValueErro('Unknown aligner `{}\''.format(name))
