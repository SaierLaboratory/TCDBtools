#!/usr/bin/env python

import argparse
import numpy as np
import sys
import sqlite3
import re
import os

try: from kevtools_common.types import TCID
except ImportError:
    TCID = None
    print('[WARNING]: Could not import kevtools_common. TCID sorting will fail if attempted.', file=sys.stderr)


def error(*things):
    print('[ERROR]:', *things, file=sys.stderr)
    exit(1)

def _parse_ranges(rangestr, index=1):
    indices = []
    for piece in rangestr.split(','):
        if piece.startswith('-'):
            if piece.count('-') >= 2: 
                startend = piece[1:].find('-')
                delim = piece[startend:].find('-')
                #endstart = piece[delim:].find('-')

                start = int(piece[:startstart])
                end = int(piece[delim+1:])
                indices.extend(list(range(start-index, end+1-index)))
            else:
                indices.append(int(piece)-index)
        elif piece.count('-') == 0:
            indices.append(int(piece)-index)
        else:
            sp = piece.split('-')
            indices.extend(list(range(int(sp[0])-index, int(sp[1])+1-index)))
        #TODO: Refactor!
    return indices

def cat(args):
    for fn in args.infile:
        with open(fn) as fh: print(fh.read())

def _scale(val, vmin, vmax):
    if val <= vmin: return 0.
    elif val >= vmax: return 1.
    else: return (val - vmin) / (vmax - vmin)

def cmap(args):
    import matplotlib.cm as cm
    colormap = cm.__getattribute__(args.cmap)
    if args.k: columns = _parse_ranges(args.k)
    else: columns = []

    if args.vmin is None or args.vmax is None:
        #bulk mode
        vmin = args.vmin
        vmax = args.vmax
        intext = ''
        for fn in args.infile:
            with open(fn) as fh:
                intext += fh.read()
        for l in intext.split('\n'):
            if l.startswith('#'): continue
            elif not l.strip(): continue
            else:
                for k, item in enumerate(l.split('\t')):
                    if (args.k is not None) and (k not in columns): continue
                    try:
                        if args.vmin is None:
                            if vmin is None: vmin = float(item)
                            else: vmin = min(vmin, float(item))
                        if args.vmax is None:
                            if vmax is None: vmax = float(item)
                            else: vmax = max(vmax, float(item))
                    except ValueError: pass
        for l in intext.split('\n'):
            if l.startswith('#'): print(l)
            elif not l.strip(): print(l)
            else:
                row = []
                for k, item in enumerate(l.split('\t')):
                    if (args.k is not None) and (k not in columns): 
                        row.append(item)
                        continue
                    try:
                        color = colormap(_scale(float(item)), vmin, vmax)[:3]
                        r = int(color[0] * 255)
                        g = int(color[1] * 255)
                        b = int(color[2] * 255)
                        row.append('\033[48;2;{};{};{}m{}\033[0m'.format(r, g, b, item))
                    except ValueError: row.append(item)
                print('\t'.join(row))

    else: 
        #line mode

        for fn in args.infile:
            with open(fn) as fh:
                for l in fh:
                    if l.startswith('#'): print(l[:-1])
                    elif not l.strip(): print(l[:-1])
                    else:
                        row = []
                        for k, item in enumerate(l[:-1].split('\t')):
                            if (args.k is not None) and (k not in columns): 
                                row.append(item)
                                continue
                            try:
                                color = colormap(_scale(float(item), args.vmin, args.vmax))[:3]
                                r = int(color[0] * 255)
                                g = int(color[1] * 255)
                                b = int(color[2] * 255)
                                row.append('\033[48;2;{};{};{}m{}\033[0m'.format(r, g, b, item))
                            except ValueError: row.append(item)
                        print('\t'.join(row))

        

def cut(args):

    neworder = []

    if args.f is None: neworder = None
    else: neworder = _parse_ranges(args.f)

    for fn in args.infile:
        fh = open(fn)
        for l in fh:
            if neworder is None:
                print(l.replace('\n', ''))
            elif not l.strip():
                print(l.replace('\n', ''))
            else:
                sl = l.replace('\n', '').split(args.d)
                fields = []
                for i in neworder:
                    if i < len(sl): fields.append(sl[i])
                    else: fields.append('')
                print('\t'.join(fields))

def hstack(args):
    if args.f is None: fill = None
    else: fill = args.f

    tables = []
    tableshapes = []

    skip = _parse_ranges(args.skip)[::-1] if args.skip else None

    #this is mostly the part where we load the raw data and compute the supertable shape
    for findex, fn in enumerate(args.infile):
        tableshapes.append([0,0])
        tables.append([])
        with open(fn) as fh:
            for l in fh:
                tableshapes[-1][0] += 1
                sl = l.replace('\n', '').split(args.d)

                if findex and skip:
                    for i in skip: sl.pop(i)

                tables[-1].append(sl)
                tableshapes[-1][1] = max(tableshapes[-1][1], len(sl))

    #if autofill is inactive, unequal table lengths is an error
    if fill is None:
        for shape in tableshapes: 
            if shape[0] != tableshapes[0][0]: 
                error('Unequal table heights')

    #else, fill in the missing values
    else:
        maxheight = max([x[0] for x in tableshapes])
        maxwidth = max([x[1] for x in tableshapes])

        for table in tables:
            for row in table:
                row.extend([fill] * (maxwidth - len(row)))

            if len(table) < maxheight:
                for r in range(maxheight - len(table)): table.append([fill] * maxwidth)

        for i in range(len(tableshapes)): tableshapes[i] = [maxheight, maxwidth]

    #now create the final table
    newtable = []
    for r in range(tableshapes[0][0]):
        newtable.append([])
        for table in tables:
            newtable[-1].extend(table[r])

    #and print to stdout
    for row in newtable: print(args.d.join(row))

def matrix(args):
    for fn in args.infile:
        if args.m: #to matrix
            matrices = {}
            k1, k2 = map(lambda x: x - 1, args.k)
            labels = []
            collabels = set()
            with open(fn) as fh:
                for lineno, line in enumerate(fh):
                    sl = line.replace('\n', '').split(args.d)
                    if args.header_line is not None:
                        if lineno < (args.header_line - 1): continue
                        elif lineno == (args.header_line - 1):
                            if line.startswith('#'): labels = line[1:].split(args.d)
                            else: labels = sl
                            continue
                    elif not labels: labels = ['k{}'.format(i) for i in range(1, len(sl)+1)]
                        
                    indices = sl[k1], sl[k2]
                    for k, val in enumerate(sl):
                        if k in (k1, k2): continue
                        elif k >= len(labels): continue #TODO: strict mode
                        field = labels[k]
                        if field not in matrices: matrices[field] = {}
                        if indices[0] not in matrices[field]: matrices[field][indices[0]] = {}
                        if indices[1] not in matrices[field][indices[0]]: matrices[field][indices[0]][indices[1]] = val
                        collabels.add(indices[1])

            for field in matrices:
                print('#{} x {} matrix of {}'.format(labels[k1], labels[k2], field))
                print('\t' + '\t'.join([str(index) for index in sorted(collabels)]))
                for row in sorted(matrices[field]):
                    print(str(row) + '\t' + '\t'.join([matrices[field][row].get(col, '') for col in sorted(collabels)]))
        elif args.t: #to columnwise table
            raise NotImplementedError

def shape(args):
    for fn in args.infile:
        
        width = 0
        height = 0
        widths = set()
        lastwidth = None
        diffs = []
        with open(fn) as fh:
            for l in fh:
                sl = l.split("\t")
                width = len(sl)
                if l.strip(): 
                    height += 1
                    widths.add(width)

                contentwidth = len([x for x in sl if len(x)])
                if lastwidth is not None: 
                    diffs.append(contentwidth - lastwidth)
                lastwidth = contentwidth

                if args.per_line:
                    print(width)

        raggedness = "[rect]"
        if len(widths) > 1: 
            if height >= 3 and (all([x >= 0 for x in diffs]) or all([x <= 0 for x in diffs])):
                raggedness = "[tri]"
            else: raggedness = "[rag]"
                
        if len(args.infile) > 1: pass
        if not args.per_line:
            print("{}{}{} {}".format(str(height).rjust(8), str(width).rjust(8), raggedness.rjust(8), fn))

def tcsort(args):
    if args.k is None: keyorder = None
    else: keyorder = _parse_ranges(args.k)

    table = []

    for fn in args.infile:
        fh = open(fn)
        for l in fh:
            if keyorder is None:
                print(l.replace('\n', ''))
            elif not l.strip():
                print(l.replace('\n', ''))
            else:
                sl = l.replace('\n', '').split(args.d)

                sortkey = []
                for k in keyorder:
                    try: sortkey.append(TCID(sl[k]))
                    except AssertionError: sortkey.append(sl[k])
                table.append((sortkey, sl))

    table.sort()
    for sk, row in table:
        print(args.d.join(row))

def transpose(args):
    table = []

    maxlen = None
    for fn in args.infile:
        fh = open(fn)
        for l in fh:
            sl = l.replace('\n', '').split(args.d)

            maxlen = max(maxlen, len(sl)) if maxlen is not None else len(sl)

            table.append(sl)

    tablet = []
    for i in range(maxlen):
        tablet.append([])
        for row in table:
            try: tablet[-1].append(row[i])
            except IndexError: tablet[-1].append('')

    if args.clone:
        for i in range(len(tablet)):
            for j in range(len(tablet[i])):
                if not len(tablet[i][j]):
                    try: tablet[i][j] = tablet[j][i]
                    except IndexError: continue

    for row in tablet:
        print(args.d.join(row))

def tsv_sql(args):
    if not args.c: 
        with open(args.o, 'w') as fh:
            fh.write('')
        return

    delimiter = '\t'
    keys = None
    rows = []

    #ONLY used for type guessing
    #no guarantee that values will retain/share indices
    cols = {}

    firstline = True
    lineno = 0
    with open(args.infile) as fh:
        for l in fh:
            lineno += 1
            if args.header_line == lineno:
                firstline = False
                if l.startswith('#'): 
                    keys = l[1:].replace('\n', '').split(delimiter)
                    continue
                else:
                    keys = l.replace('\n', '').split(delimiter)
                    continue
            elif firstline:
                firstline = False
                if l.startswith('#'): 
                    keys = l[1:].replace('\n', '').split(delimiter)
                    continue
            if not l.strip(): continue
            else:
                if l.startswith('#'): continue
                row = l.replace('\n', '').split(delimiter)
                for i, item in enumerate(row):
                    if i not in cols: 
                        cols[i] = []
                        if rows:
                            for j in range(len(rows)):
                                rows[j] = rows[j] + (None,) * (len(cols) - len(rows[j]))
                    cols[i].append(item)
                if len(row) < len(cols):
                    if args.ignore_errors:
                        row.extend([None] * (len(cols) - len(row)))
                    else: error('Incomplete row at line {}. (Use --ignore-errors to suppress)'.format(lineno))
                rows.append(tuple(row))

    if args.ignore_column_names or (keys is None):
        keys = ['c{}'.format(i+1) for i in cols]

    types = {}
    for i in cols:
        lasttype = None
        types[i] = str
        for x in cols[i][:100]:
            if lasttype is None: lasttype = _guess_type(x)
            elif lasttype is str: 
                types[i] = str
                break
            elif lasttype != _guess_type(x): 
                types[i] = str
                break
        types[i] = lasttype

    sanitized_fn = os.path.basename(args.infile)
    if '.' in sanitized_fn: 
        sanitized_fn = '.'.join(os.path.splitext(sanitized_fn)[:-1])#sanitized_fn[:sanitized_fn.find('.')]
    #sanitized_fn = 'table_' + re.sub('[^-_a-zA-Z0-9]', '_', sanitized_fn)
    sanitized_fn = 'tsv'

    db = sqlite3.connect(':memory:')
    curs = db.cursor()
    columnstr = ''
    typenames = {int:'INTEGER', float:'FLOAT', str:''}
    for key, typeindex in zip(keys, types):
        dtype = types[typeindex]
        columnstr += '{} {}, '.format(key, typenames[dtype])
    if columnstr.endswith(', '): columnstr = columnstr[:-2]
    createcmd = "CREATE TABLE {tablename}({columnstr});".format(tablename=sanitized_fn, columnstr=columnstr)
    curs.execute(createcmd)

    curs.executemany("INSERT INTO {tablename}({keys}) VALUES ({questionmarks})".format(tablename=sanitized_fn, keys=', '.join(keys), questionmarks=','.join(['?' for x in keys])), rows)

    try: curs.execute(args.c)
    except sqlite3.OperationalError as e:
        error('SQL error: {}'.format(e))

    with open(args.o, 'w') as fh:
        for row in curs.fetchall():
            fh.write(delimiter.join([str(x) for x in row]) + '\n')
        

def _guess_type(s):
    try: 
        int(s)
        return int
    except ValueError:
        try:
            float(s)
            return float
        except ValueError: 
            if '%' in s:
                float(s[:-1])
                return float
            elif '$' in s:
                float(s[1:])
                return float
            else: return str
            return str
    return str


def uniq(args):
    if args.k is None: keyorder = None
    else: keyorder = _parse_ranges(args.k)

    done = set()
    
    for fn in args.infile:
        with open(fn) as fh:
            for l in fh:
                sl = l.replace('\n', '').split(args.d)

                if keyorder is None: relevant = l
                else: relevant = tuple([sl[i] for i in keyorder])

                if relevant in done: pass
                else: 
                    print(l.replace('\n', ''))
                    done.add(relevant)

###

def get_full_parser():
    parser = argparse.ArgumentParser()

    #parser.add_argument('command', nargs='*', default='cat', help='What to do with files. Enter help for details')

    subparsers = parser.add_subparsers(help='tsviss army knife commands')


    ####cat###
    parser_cat = subparsers.add_parser('cat', help='Concatenate files. Utterly useless command that should never have been implemented.')
    parser_cat.set_defaults(func=cat)

    parser_cat.add_argument('infile', 
        nargs='*', default=['/dev/stdin'], 
        help='Files to process')

    #parser_cat.add_argument('infile', 
    #    nargs='*', default=['/dev/stdin'], 
    #    help='Files to concatenate. Defaults to stdin.')

    #parser_cat.add_argument('-b',
    #    action='store_true',
    #    help='Number non-blank output lines. 1-indexed.')

    #parser_cat.add_argument('-n',
    #    action='store_true',
    #    help='Number output lines. 1-indexed.')

    #parser_cat.add_argument('--no-comment',
    #    action='store_true',
    #    help='Remove all #-initial lines')

    #parser_cat.add_argument('--no-empty',
    #    action='store_true',
    #    help='Remove all empty/whitespace-only lines')

    ###color###
    parser_cmap = subparsers.add_parser('cmap', help='Apply conditional coloring to TSVs with numeric data. Requires matplotlib.')
    parser_cmap.set_defaults(func=cmap)

    parser_cmap.add_argument('infile', 
        nargs='*', default=['/dev/stdin'], 
        help='Files to process')

    parser_cmap.add_argument('--cmap',
        default='viridis',
        help='Color map to use. Default: viridis')

    parser_cmap.add_argument('-k',
        help='Comma-separated list of column ranges to color. 1-indexed.')

    parser_cmap.add_argument('--vmin',
        type=float,
        help='Manually set the "low" end of the gradient. Switches to line-by-line mode if used in conjunction with --vmax which may be desirable if your dataset is large.')

    parser_cmap.add_argument('--vmax',
        type=float,
        help='Manually set the "high" end of the gradient. Switches to line-by-line mode if used in conjunction with --vmin which may be desirable if your dataset is large.')

    ###cut###
    parser_cut = subparsers.add_parser('cut', help='Cut files. Allows duplicate/non-increasing column IDs.')
    parser_cut.set_defaults(func=cut)

    parser_cut.add_argument('-d', 
        default='\t',
        help='Field delimiter. Defaults to <TAB>')

    parser_cut.add_argument('-f',
        default=None,
        help='Which fields to keep')

    parser_cut.add_argument('-s',
        action='store_true',
        help='Remove undelimited lines instead')

    parser_cut.add_argument('infile', 
        nargs='*', default=['/dev/stdin'], 
        help='Files to process')

    ###hstack###
    parser_hstack = subparsers.add_parser('hstack', help='Ligate TSVs horizontally. Crashes by default if table heights differ or if non-rightmost tables have ragged right edges.')
    parser_hstack.set_defaults(func=hstack)

    parser_hstack.add_argument('-d', 
        default='\t',
        help='Field delimiter. Defaults to <TAB>')

    parser_hstack.add_argument('-f',
        default=None,
        help='Fill empty cells with this character if set. Otherwise, mismatched table heights will result in crashes.')

    parser_hstack.add_argument('--skip',
        help='Columns to exclude in non-initial tables. Useful for skipping identical row labels provided all tables are in the same order.')

    parser_hstack.add_argument('infile', 
        nargs='*', default=['/dev/stdin'], 
        help='Files to process')

    ###matrix###
    parser_matrix = subparsers.add_parser('matrix', help='Convert TSVs into matrices or vice versa.')
    parser_matrix.set_defaults(func=matrix)

    parser_matrix.add_argument('-d',
        default='\t',
        help='Field delimiter. Defaults to <TAB>')

    parser_matrix.add_argument('--header-line', 
        type=int, default=None, 
        help='Line to use as a header. If set, tsviss_army_knife will skip lines preceding the header line.')

    parser_matrix_direction = parser_matrix.add_mutually_exclusive_group(required=True)
    parser_matrix_direction.add_argument('-m',
        action='store_true',
        help='Convert a TSV into a matrix or series of matrices')

    parser_matrix.add_argument('-k',
        type=int, nargs=2, default=[1, 2],
        help='1-indexed columns to use as matrix indices in row-column order (default: 1, 2)')

    parser_matrix_direction.add_argument('-t',
        action='store_true',
        help='Convert a matrix or series of matrices into a TSV')

    parser_matrix.add_argument('infile', 
        nargs='*', default=['/dev/stdin'], 
        help='Files to process')

    ###shape###
    parser_shape = subparsers.add_parser('shape', help='Get TSV shape stats.')
    parser_shape.set_defaults(func=shape)

    parser_shape.add_argument('--per-line',
        action='store_true',
        help='Output per-line stats')
    parser_shape.add_argument('infile', 
        nargs='*', default=['/dev/stdin'], 
        help='Files to process')

    ###sql###
    parser_sql = subparsers.add_parser('sql', help='Run SQLite commands against TSVs')
    parser_sql.add_argument('infile', 
        nargs='?',
        default='/dev/stdin',
        help='File to process (default:stdin)')

    parser_sql.add_argument('-c',
        help='SQLite commands to execute. The table of imported data is named "tsv". If the first line in the input isn\'t a #-commented row of column names, tsviss_army_knife will name the columns c1, c2, c3, ... cn.')

    parser_sql.add_argument('-o',
        default='/dev/stdout',
        help='Output file (default: stdout)')

    parser_sql.add_argument('--ignore-errors',
        action='store_true',
        help='Ignore row size mismatch errors')

    parser_sql.add_argument('--header-line', type=int, default=None, help='Line to use as a header. If set, tsviss_army_knife will skip lines preceding the header line.')

    parser_sql.add_argument('--ignore-column-names', action='store_true', help='Ignore column names in header if available in favor of c1, c2, c3...')

    parser_sql.set_defaults(func=tsv_sql)

    ###tcsort###
    parser_tcsort = subparsers.add_parser('tcsort', help='Sort files by TC-ID. Non-TCIDs are sorted alphabetically and lower than TC-IDs')
    parser_tcsort.set_defaults(func=tcsort)

    parser_tcsort.add_argument('-d', 
        default='\t',
        help='Field delimiter. Defaults to <TAB>')

    parser_tcsort.add_argument('-k', 
        default=None,
        help='Comma-separated list of ranges of 1-indexed columns with TC-IDs')

    parser_tcsort.add_argument('infile', 
        nargs='*', default=['/dev/stdin'], 
        help='Files to process')

    ###transpose###
    parser_transpose = subparsers.add_parser('transpose', help='Transpose tables.')
    parser_transpose.set_defaults(func=transpose)

    parser_transpose.add_argument('-d', 
        default='\t',
        help='Field delimiter. Defaults to <TAB>')

    parser_transpose.add_argument('--clone',
        action='store_true',
        help='Fill empty fields from transpose')

    parser_transpose.add_argument('infile', nargs='*', default=['/dev/stdin'], help='Files to process')

    ###uniq###
    parser_uniq = subparsers.add_parser('uniq', help='Collect unique entries. Unlike POSIX uniq, this picks up on multicolumn uniqueness keys (e.g. qacc and sacc in BLAST TSVs) and nonconsecutive duplicates.')
    parser_uniq.set_defaults(func=uniq)

    parser_uniq.add_argument('-c', 
        default='\t',
        help='Count occurrences')

    parser_uniq.add_argument('-d', 
        default='\t',
        help='Field delimiter. Defaults to <TAB>')

    parser_uniq.add_argument('-k', 
        default=None,
        help='Comma-separated list of ranges of 1-indexed columns with TC-IDs')

    parser_uniq.add_argument('infile', 
        nargs='*', default=['/dev/stdin'], 
        help='Files to process')

    return parser

if __name__ == '__main__':
    parser = get_full_parser()
    args = parser.parse_args()

    if getattr(args, "func", None) is None:
        parser.print_help()
        exit(1)
    else: args.func(args)
