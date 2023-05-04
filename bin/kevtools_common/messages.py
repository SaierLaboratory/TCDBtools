import sys

VERBOSITY = 3
#0: No messages
#1: Errors only
#2: Errors and warnings
#3: All messages

def info(*things):
	if VERBOSITY >= 3: print('[INFO]', *things, file=sys.stderr)

def warn(*things):
	if VERBOSITY >= 2: print('[WARNING]', *things, file=sys.stderr)

def error(*things):
	if VERBOSITY >= 1: print('[ERROR]', *things, file=sys.stderr)
	exit(1)
