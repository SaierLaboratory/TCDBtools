import re

def parse_hmmtop(top):
    out = []
    if isinstance(top, bytes): top = top.decode('ascii')
    for l in top.split('\n'):
        if not l.strip(): continue
        indices = re.findall('(?:[0-9]+ *)+$', l)[0]
        indices = [int(x) for x in indices.strip().split()[1:]]
        out.append([[indices[i]-1, indices[i+1]] for i in range(0, len(indices), 2)])
    return out

class Mapper(object):
    ''' Compact index mapper for sequences '''
    def __init__(self): 
        self.offsets = {}
        self._sorted_indices = []

    def __getitem__(self, index):
        for k in self._sorted_indices:
            if k <= index: index += self.offsets[k]
            elif k > index: break
        return index

    def get(self, index): return self.__getitem__(index)

    def __setitem__(self, index, value):
        self.offsets[index] = value
        self._sorted_indices.append(index)
        self._sorted_indices.sort()

    def __repr__(self):
        return 'Mapper({})'.format(self.offsets)

def trueresi_to_hmmtop(seq):
    mapping = Mapper()
    for resi, resn in enumerate(seq):
        if resn not in 'ACDEFGHIKLMNPQRSTVWY': mapping[resi] -= 1
    return mapping

def hmmtop_to_trueresi(seq):
    mapping = Mapper()
    for resi, resn in enumerate(seq):
        if resn not in 'ACDEFGHIKLMNPQRSTVWY': mapping[resi] += 1
    return mapping

def parse_ranges(s, spandelimiter=',', endpointdelimiter='-', strict=True, open=False, ranges=False, offset=0):
    '''
    Arguments:

    s: String to parse into ranges

    spandelimiter: Substring to use as a delimiter between ranges
    endpointdelimiter: Substring to use as a delimiter between endpoints
    strict: Raise ValueErrors for integer endpoints and null intervals that end before they start
    open: Output open intervals, i.e. intervals containing their endpoints, instead of half-open (left-open, right-closed) intervals. This is not Pythonic even if it makes sense when producing human-readable output which is why it is not enabled by default.
    ranges: Output a list of ranges instead of a list of 2-element lists. Ignores the open argument in order to generate ranges that contain the given endpoints
    offset: Offset to apply to all indices
    '''
    spans = []
    for interval in s.split(spandelimiter):
        values = interval.split(endpointdelimiter)

        if not values: continue

        if strict: endpoints = [int(x) for x in values]
        else: 
            for x in values:
                try: endpoints.append(int(x))
                except ValueError: warnings.warn('Interval {} in string {} contains an invalid endpoint, skipping'.format(interval, s))

        if len(values) == 1:
            endpoints = endpoints * 2
        elif strict and len(endpoints) != 2:
            raise ValueError('Incorrect number of endpoints for range {} in string {}'.format(interval, s))

        #for half-open residues, stretch the interval so it covers the right endpoint per the convention that "1-1" ::= [1, 1] and not [1, 1)
        #as these are integer intervals, [1, 2) is equivalent to [1, 1]
        if ranges or not open: endpoints[1] += 1

        if offset: endpoints = list(map(lambda x: x + offset, endpoints))

        if ranges: spans.append(range(*endpoints))
        else: spans.append(endpoints)

    return spans


