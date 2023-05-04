#!/usr/bin/env python

import argparse 
import configparser
import shlex
import libquod
import re
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.transforms as transforms
import matplotlib
import os
import sys
from Bio import SeqIO, AlignIO

STRICT = False

IMPLICIT_ENTITIES = set((
    'title', 'ltitle', 'ctitle', 'rtitle',
    'label', 'xlabel', 'ylabel',
    'ticks', 'xticks', 'yticks'
))
def warn(*things): 
    print('[WARNING]', *things, file=sys.stderr)
def error(*things): 
    print('[ERROR]', *things, file=sys.stderr)
    exit(1)

class SetOfAllThings(set):
    #sorry, Russell
    def __contains__(self, other): return True

class _Parsing(object):
    allstrs = frozenset(['all', 'yes', 'true'])
    nonestrs = frozenset(['no', 'none', 'false'])

    @staticmethod
    def get_list(cfgstr, dtype=None): 
        if dtype is None: return shlex.split(cfgstr)
        else: return [dtype(x) for x in shlex.split(cfgstr)]

    @staticmethod
    def get_seqlist(cfgstr):
        filenames = []
        flags = []
        #0: full sequence
        #1: next sequence is a fragment
        #2: next sequence is a full sequence (hidden)
        #3: next sequence is a full sequence
        flag = 0
        nextflag = [0, 2, 3, 0]
        for token in shlex.split(cfgstr):
            flag = nextflag[flag]
            if token == "--frag":
                flag = 1
            else:
                filenames.append(token)
                flags.append(flag - 1 if flag else flag)

        return filenames, flags
        

    @staticmethod
    def get_listlist(cfgstr): return [shlex.split(line) for line in cfgstr.split('\n')]

    @staticmethod
    def apply_selector(cfgstr, arr, one_indexed=False):
        if cfgstr.lower() in _Parsing.allstrs: return arr[:]
        elif cfgstr.lower() == _Parsingnonestrs: return arr[-1:0]
        #TODO: Find out whether always outputting a list would be a good idea for the above

        elif '-' in cfgstr: 
            ranges = _Parsing.get_list(cfgstr)
            indices = []
            for r in ranges:
                if '-' in r: 
                    start, end = r.split('-')
                    start = int(start)
                    end = int(end)

                    if one_indexed: start -= 1
                    else: end += 1
                    indices.extend(range(start, end))
            return [arr[i] for i in indices]
        else:
            if one_indexed: return [arr[int(i) - 1] for i in _Parsing.get_list(cfgstr)]
            else: return [arr[int(i)] for i in _Parsing.get_list(cfgstr)]

    @staticmethod
    def load_sequences(fn):
        if fn.endswith('.clu') or fn.endswith('.aln'): return AlignIO.parse(fn, 'clustal')

        else: return SeqIO.parse(fn, 'fasta')

    @staticmethod
    def load_msa(fn):
        pass

    @staticmethod
    def get_regions(cfgstr):
        entlist = []

        for l in cfgstr.split('\n'):
            if not l.strip(): continue
            elif l.strip().startswith('#'): continue

            else:
                obj = libquod.entities.Region()

                tokens = shlex.split(l)
                for t in tokens:
                    if t.startswith('label:') or t.startswith('text:'): 
                        obj.text = t[t.find(':')+1:]
                    elif t.startswith('fontsize:'): 
                        obj.fontsize = int(t[t.find(':')+1:])
                    elif t.startswith('ha'):
                        obj.halign = t[t.find(':')+1:]
                    elif t.startswith('va'):
                        obj.valign = t[t.find(':')+1:]
                    elif t.startswith('textcolor'):
                        obj.textcolor = t[t.find(':')+1:]
                    elif t.startswith('alpha'):
                        obj.alpha = float(t[t.find(':')+1:])
                    elif t.startswith('color') or t.startswith('fc') or t.startswith('facecolor'):
                        obj.facecolor = t[t.find(':')+1:]
                    elif t.startswith('x:'):
                        spanstr = t[t.find(':')+1:]
                        obj.spans = [[float(x) for x in span.split(':')] for span in spanstr.split(',')]
                    elif t.startswith('y:'):
                        obj.yspan = [float(y) for y in t[t.find(':')+1:].split(':')]
                    elif t.startswith('period:'):
                        obj.period = float(t[t.find(':')+1:])
                    elif t.startswith('dashed:'):
                        obj.dashed = float(t[t.find(':')+1:])
                    elif t.startswith('offset:'):
                        obj.offset = float(t[t.find(':')+1:])
                    else: error('Unrecognized token in region: "{}"'.format(t))
                entlist.append(obj)

        return entlist

    @staticmethod
    def get_walls(cfgstr):
        entlist = []

        for l in cfgstr.split('\n'):
            if not l.strip(): continue
            elif l.strip().startswith('#'): continue

            else:
                obj = libquod.entities.Wall()

                tokens = shlex.split(l)
                for t in tokens:
                    if t.startswith('label:') or t.startswith('text:'): 
                        obj.text = t[t.find(':')+1:]
                    elif t.startswith('fontsize:'): 
                        obj.fontsize = int(t[t.find(':')+1:])
                    elif t.startswith('ha'):
                        obj.halign = t[t.find(':')+1:]
                    elif t.startswith('va'):
                        obj.valign = t[t.find(':')+1:]
                    elif t.startswith('textcolor'):
                        obj.textcolor = t[t.find(':')+1:]
                    elif t.startswith('alpha'):
                        obj.alpha = float(t[t.find(':')+1:])
                    elif t.startswith('color') or t.startswith('fc') or t.startswith('facecolor'):
                        obj.facecolor = t[t.find(':')+1:]
                    elif t.startswith('x:'):
                        spanstr = t[t.find(':')+1:]
                        #obj.spans = [[float(x) for x in span.split(':')] for span in spanstr.split(',')]
                        obj.spans = libquod.entities.Parser.parse_ranges(spanstr, allow_points=True)
                    elif t.startswith('y:'):
                        obj.yspan = [float(y) for y in t[t.find(':')+1:].split(':')]
                        obj.y = float(t[t.find(':')+1:])
                    elif t.startswith('ypos:'):
                        obj.ypos = t[t.find(':')+1:]
                    elif t.startswith('scale:'):
                        obj.scale = float(t[t.find(':')+1:])
                    else: error('Unrecognized token in region: "{}"'.format(t))
                if obj.y is None:
                    if obj.ypos == '-': obj.y = -1.732
                    else: obj.y = 1.732
                entlist.append(obj)

        return entlist

    @staticmethod
    def get_tmss(cfgstr):
        entlist = []

        for l in cfgstr.split('\n'):
            if not l.strip(): continue
            elif l.strip().startswith('#'): continue

            else:
                obj = libquod.entities.HMMTOP()
                savetext = False
                textkwargs = {}

                tokens = shlex.split(l)
                for t in tokens:
                    if t.startswith('label:') or t.startswith('text:'): 
                        textkwargs['text'] = t[t.find(':')+1:]
                        savetext = True
                    elif t.startswith('fontsize:'): 
                        textkwargs['fontsize'] = int(t[t.find(':')+1:])
                    elif t.startswith('ha:'):
                        textkwargs['halign'] = t[t.find(':')+1:]
                    elif t.startswith('va:'):
                        textkwargs['valign'] = t[t.find(':')+1:]
                    elif t.startswith('textcolor'):
                        textkwargs['textcolor'] = t[t.find(':')+1:]
                    elif t.startswith('alpha'):
                        obj.alpha = float(t[t.find(':')+1:])
                    elif t.startswith('color:') or t.startswith('fc:') or t.startswith('facecolor:'):
                        obj.facecolor = t[t.find(':')+1:]
                    elif t.startswith('x:'):
                        spanstr = t[t.find(':')+1:]
                        obj.spans = [[float(x) for x in span.split(':')] for span in spanstr.split(',')]
                    elif t.startswith('y:'):
                        textkwargs['pos'][1] = float(t[t.find(':')+1:])
                    else: error('Unrecognized token in tms: "{}"'.format(t))
                entlist.append(obj)
                if savetext: pass

        return entlist

class MultiquodConfig(object):
    dpi = None
    equal = True
    tight = True
    width = None
    height = None
    outfn = None
    props = None

    def __init__(self):
        self.subplots = {} #? or maybe list?
        self.defaults = None
        self.entities = {}
        self.rows = []
        self.rowsubplots = []

        self.props = {}

        self.config = configparser.ConfigParser()
        self.initialize_config()

    #useful mutators

    def add_row(self, row=None):
        if row is None: self.rows.append([])
        else: self.rows.append(row)


    def add_subplot(self, subplot=None):
        if subplot is None: 
            name = 'subplot{}'.format(len(self.subplots))
            subplots[name] = Subplot(name=name)


    def render(self):
        for name in self.subplots: self.subplots[name].render()

        self._allocate_grid()


    def _allocate_grid(self):
        ''' An evil method whose only purpose is to produce arbitrarily sized plots '''

        for row in self.rows:
            lengths = []
            for name in row: 
                length = None
                if self.subplots[name].length: lengths.append(self.subplots[name].length)
                else: 
                    try: lengths.append(self.subplots[name].xlim[1] - self.subplots[name].xlim[0])
                    except TypeError: lengths.append(20)

    #IO methods

    def write(self, fh):
        out = ''
        fh.write(out.encode('utf-8'))


    @staticmethod
    def load(fn):
        if not os.path.isfile(fn): raise FileNotFoundError('Could not find file `{}\''.format(fn))
        obj = MultiquodConfig()

        obj.config.read(fn)

        obj.process_config()

        return obj


    def initialize_config(self):
        self.config['figure'] = {}
        self.config['subplots.default'] = {}

        #self.config['subplots.default'] = dict(
        #            ctitle='',
        #            legend='no',
        #            ltitle='',
        #            rtitle='',
        #            title='',
        #            xlabel='',
        #            xlim='',
        #            xscale=1,
        #            ylabel=None,
        #            ylim=None,
        #            yscale=1,
        #            length=None,
        #            mode='',
        #            kernel='flat',
        #        )


    def draw(self, fig=None, quiet=False):
        if fig is None: fig = plt.figure()

        if self.height: fig.set_figheight(self.height, forward=True)
        if self.width: fig.set_figwidth(self.width, forward=True)

        if self.tight: fig.set_tight_layout(True)


        if len(self.subplots) == 1:
            for name in self.subplots: subplot = self.subplots[name]
            print(id(fig))
            subplot.ax = fig.add_subplot(111)
            self.rowsubplots = [subplot.ax]
            
            subplot.draw()
            subplot.ax.set_xlim(subplot.xlim)
            subplot.ax.set_ylim(subplot.ylim)

            subplot.ax.axhline(0, color='k', lw=1)

        elif len(self.subplots) == len(self.rows): 
            gs = gridspec.GridSpec(len(self.rows), 1, figure=fig)

            for rowi, row in enumerate(self.rows):
                subplot = self.subplots[row[0]]
                #self.rowsubplots.append(fig.add_subplot(
                subplot.ax = fig.add_subplot(gs[rowi,:])
            
                subplot.draw()
                subplot.ax.set_xlim(subplot.xlim)
                subplot.ax.set_ylim(subplot.ylim)
                subplot.ax.axhline(0, color='k', lw=1)

        elif self.equal: #aka hvordanlike
            grids = {}
            for row in self.rows: 
                if len(row) in grids: continue
                grids[len(row)] = gridspec.GridSpec(len(self.rows), len(row), figure=fig)

            for rowi, row in enumerate(self.rows):
                gs = grids[len(row)]
                for subploti, name in enumerate(row):
                    subplot = self.subplots[name]
                    subplot.ax = fig.add_subplot(gs[rowi,subploti])

                    subplot.draw()
                    subplot.ax.set_xlim(subplot.xlim)
                    subplot.ax.set_ylim(subplot.ylim)
                    subplot.ax.axhline(0, color='k', lw=1)

        else: #:( what a nightmare!
            gs = gridspec.GridSpec(len(self.rows), 1, figure=fig)
            if self.tight: gs.tight_layout(fig)
            for rowi, row in enumerate(self.rows):
                xlimlist = []
                for subploti, name in enumerate(row):
                    subplot = self.subplots[name]
                    if subplot.length is not None: xlimlist.append(subplot.length)
                    else:
                        xlim, ylim = subplot.get_tightbbox()

                        if xlim[0] is None or xlim[1] is None: xlimlist.append(None)
                        else: xlimlist.append(xlim[1] - xlim[0])

                if None in xlimlist: raise ValueError

                #FIXME: this probably breaks horribly on metric systems!
                rowspec = gridspec.GridSpecFromSubplotSpec(1, len(row), subplot_spec=gs[rowi], width_ratios=xlimlist, wspace=45/72/6.4**2*fig.get_figwidth())

                for subploti, name in enumerate(row):
                    subplot = self.subplots[name]
                    subplot.ax = fig.add_subplot(rowspec[:, subploti:subploti+1])

                    subplot.draw()
                    subplot.ax.set_xlim(subplot.xlim)
                    subplot.ax.set_ylim(subplot.ylim)
                    subplot.ax.axhline(0, color='k', lw=1)


        for name in self.subplots: 
            subplot = self.subplots[name]
            if subplot.xscale is not None:
                xlim = subplot.ax.get_xlim()

                if xlim[1] % subplot.xscale: xmax = xlim[1]
                else: xmax = xlim[1] + subplot.xscale

                subplot.ax.set_xticks(np.arange(xlim[0], xmax, subplot.xscale))

            if subplot.yscale is not None:
                ylim = subplot.ax.get_ylim()

                if ylim[1] % subplot.yscale: ymax = ylim[1]
                else: ymax = ylim[1] + subplot.yscale

                subplot.ax.set_yticks(np.arange(ylim[0], ylim[1]+subplot.yscale, subplot.yscale))

            else:
                
                ylim = subplot.ax.get_ylim()

                if ylim[1] % 1: ymax = ylim[1]
                else: ymax = ylim[1] + 1

                subplot.ax.set_yticks(np.arange(ylim[0], ylim[1]+1, 1))

            if subplot.legend: subplot.ax.legend()


        if self.outfn is not None: fig.savefig(self.outfn, dpi=self.dpi)


        #FIXME: don't always show
        if not quiet: plt.show()

        #FIXME: don't always close in case the user wants to make further modifications to fig
        else: 
            plt.close()

        return fig
        

    #parsing methods

    def process_config(self):
        hasdefault = False
        #preprocessing step
        for section in self.config.sections():
            if section.startswith('subplots.'): 
                name = section[section.find('.')+1:]
                if name.startswith('default'):
                    hasdefault = True
                    self.defaults = Subplot(name=name)
                    self.defaults.process_config(self.config[section])
                else:
                    self.subplots[name] = Subplot(name=name)


        name = section[section.find('.')+1:]
        if hasdefault and name != 'default':
            #FIXME: find a less horrible way to do this
            #excluded for simplicity: entities, mplobjects, externs, externflags

            for attrib in ('modes', 'props', 'ltitle', 'ctitle', 'rtitle', 'title', 'legend', 'xlim', 'ylim', 'xlabel', 'ylabel', 'length', 'notms'):

                value = self.defaults.__getattribute__(attrib)
                try: 
                    hash(value)
                    hashable = True
                except TypeError: hashable = False

                if type(value) is list:
                    setattr(self.subplots[name], attrib, self.defaults.__getattribute__(attrib)[:])
                elif type(value) is dict:
                    newdict = {}
                    newdict.update(self.defaults.__getattribute__(attrib))
                    setattr(self.subplots[name], attrib, newdict)
                elif hashable:
                    setattr(self.subplots[name], attrib, self.defaults.__getattribute__(attrib))
                else:
                    raise TypeError('Undefined copy method for default attribute "{}"'.format(attrib))

            #self.subplots[name].process_config(self.config[section])


        for section in self.config.sections():
            if section == 'figure':
                for key in self.config[section]:
                    if key == 'rows': 
                        self.rows = _Parsing.get_listlist(self.config[section][key])
                        for row in self.rows:
                            for name in row:
                                if name not in self.subplots: self.subplots.update({name:Subplot(name)})
                    elif key == 'dpi':
                        self.dpi = self.config[section].getint('dpi')
                    elif key == 'equal':
                        self.equal = self.config[section].getboolean('equal')
                    elif key == 'height':
                        self.height = self.config[section].getfloat('height')
                    elif key == 'save':
                        self.outfn = self.config[section].get('save')
                    elif key == 'tight':
                        self.tight = self.config[section].getboolean('tight')
                    elif key == 'width':
                        self.width = self.config[section].getfloat('width')

                    elif '.' in key: 
                        entname, propname = key.split('.')
                        if propname == 'font': self.add_prop(entname, propname, self.config[section].getfloat(key))
                        elif propname == 'color': self.add_prop(entname, propname, self.config[section].get(key))
                        elif propname == 'linewidth': self.add_prop(entname, propname, self.config[section].getfloat(key))
                        elif 'size' in propname: self.add_prop(entname, propname, self.config[section].getfloat(key))
                        elif 'alpha' in propname: self.add_prop(entname, propname, self.config[section].getfloat(key))
                        else: self.add_prop(entname, propname, self.config[section].get(key))
                        
                    elif not STRICT: warn('Unrecognized key in [layout]: "{}"'.format(key))

                    else: error('Unrecognized key in [layout]: "{}"'.format(key))

            elif section.startswith('subplots.'):
                name = section[section.find('.')+1:]
                if not name.startswith('default'):
                    self.subplots[name].process_config(self.config[section])
                else: 
                    self.defaults = Subplot(name='default')
                    self.defaults.process_config(self.config[section])

        if self.subplots and not self.rows:
            for name in self.subplots: self.rows.append([name])

        #apply global props where necessary
        for entname in self.props:
            for propname in self.props[entname]:
                value = self.props[entname][propname]

                for subplotname in self.subplots:
                    subplot = self.subplots[subplotname]
                    subplot.add_prop(entname, propname, value, force=False)

    def add_prop(self, entname, propname, value):
        if entname not in self.props: self.props[entname] = {}
        if propname not in self.props[entname]: self.props[entname][propname] = {}
        self.props[entname][propname] = value


    #comparison methods

    def __eq__(self, other):
        if not isinstance(other, type(self)): return False
        else:
            try: 
                return (self.subplots == other.subplots) and (self.entities == other.entities) and (self.rows == other.rows)
            except AttributeError: return False

class Subplot(object):
    def __init__(self, name=None, entities=None, ax=None):
        self.name = name
        self.entities = {} if entities is None else entities
        self.mplobjects = {}
        self.externs = []
        self.externflags = []
        self.modes = []
        self.ax = ax

        self.props = {}

        self.ltitle = None
        self.ctitle = None
        self.rtitle = None
        self.title = None

        self.legend = False

        self.xlim = None
        self._lock_xlim = False
        self.xscale = None
        self.ylim = None
        self._lock_ylim = False
        self.yscale = None

        self.xlabel = None
        self.ylabel = None

        self.length = None
        self._lock_length = False

        self.notms = []

    def __repr__(self):
        return '{type}(name={name})'.format(type=type(self).__name__, name=self.name)

    def dump(self):
        obj = {}
        for attrib in dir(self):
            obj[attrib] = self.__getattribute__(attrib)
        return obj

    def get_tightbbox(self):
        totalxlim = [None, None]
        totalylim = [None, None]

        for entname in self.entities:
            entity = self.entities[entname]
            xlim, ylim = entity.get_bounding_box()

            #FIXME: unify all lim definitions
            if xlim is None: pass
            elif totalxlim[0] is None: totalxlim[0] = xlim[0]
            elif xlim[0] is None: pass
            else: totalxlim[0] = min(totalxlim[0], xlim[0])

            if ylim is None: pass
            elif totalylim[0] is None: totalylim[0] = ylim[0]
            elif ylim[0] is None: pass
            else: totalylim[0] = min(totalylim[0], ylim[0])

            if xlim is None: pass
            elif totalxlim[1] is None: totalxlim[1] = xlim[1]
            elif xlim[1] is None: pass
            else: totalxlim[1] = min(totalxlim[1], xlim[1])

            if ylim is None: pass
            elif totalylim[1] is None: totalylim[1] = ylim[1]
            elif ylim[1] is None: pass
            else: totalylim[1] = max(totalylim[1], ylim[1])

        #FIXME: This will cause problems with stuff down the line! Use a sane default like [0,1], [0,1] next time!
        if self._lock_xlim: totalxlim = self.xlim
        if self._lock_ylim: totalylim = self.ylim

        return totalxlim, totalylim
            

    def draw(self, ax=None):

        if ax is None: 
            if self.ax is None: raise TypeError('No ax defined for subplot')
            else: ax = self.ax
        else: ax = ax
        self.ax = ax

        for entname in self.entities:
            self.mplobjects[entname] = self.entities[entname].plot(ax)

        self.apply_title_props()
            
        self.apply_label_props()
            
        self.apply_tick_props()

    def apply_label_props(self):
        xparams = {}
        yparams = {}
        xparams.update(self.props.get('label', {}))
        xparams.update(self.props.get('xlabel', {}))
        yparams.update(self.props.get('label', {}))
        yparams.update(self.props.get('ylabel', {}))

        if self.xlabel is not None: self.ax.set_xlabel(self.xlabel.replace('<SPACE>', ' '), xparams)
        if self.ylabel is not None: self.ax.set_ylabel(self.ylabel.replace('<SPACE>', ' '), yparams)

    def apply_tick_props(self):
        self.ax.tick_params(axis='both', **self.props.get('ticks', {}))
        self.ax.tick_params(axis='x', **self.props.get('xticks', {}))
        self.ax.tick_params(axis='y', **self.props.get('yticks', {}))

    def apply_title_props(self):
        if self.title is None and self.ltitle is None and self.ctitle is None and self.rtitle is None: return
        if self.ctitle is None and self.title is not None: self.ctitle = self.title

        titlepad = matplotlib.rcParams['axes.titlepad']
        transOffset = transforms.ScaledTranslation(0, titlepad/72, self.ax.figure.dpi_scale_trans)

        titlefont = ltitlefont = ctitlefont = rtitlefont = None

        if 'title' in self.props and 'fontsize' in self.props['title']: 
            titlefont = ltitlefont = ctitlefont = rtitlefont = self.props['title']['fontsize']

        if 'ltitle' in self.props and 'fontsize' in self.props['ltitle']: 
            ltitlefont = self.props['ltitle']['fontsize']
        if 'ctitle' in self.props and 'fontsize' in self.props['ctitle']: 
            ctitlefont = self.props['ctitle']['fontsize']
        if 'rtitle' in self.props and 'fontsize' in self.props['rtitle']: 
            rtitlefont = self.props['rtitle']['fontsize']

        if self.ltitle is not None:
            t = self.ax.text(0.0, 1.0, self.ltitle,
                ha='left', va='baseline', transform=self.ax.transAxes, fontsize=ltitlefont)
            t.set_transform(self.ax.transAxes + transOffset)
        if self.ctitle is not None:
            t = self.ax.text(0.5, 1.0, self.ctitle,
                ha='center', va='baseline', transform=self.ax.transAxes, fontsize=ctitlefont)
            t.set_transform(self.ax.transAxes + transOffset)
        if self.rtitle is not None: 
            t = self.ax.text(1.0, 1.0, self.rtitle,
                ha='right', va='baseline', transform=self.ax.transAxes, fontsize=rtitlefont)
            t.set_transform(self.ax.transAxes + transOffset)


    def process_config(self, cfgdict=None, defaultcfgdict=None):
        if cfgdict is None: return

        if defaultcfgdict is not None: 
            for key in defaultcfgdict:
                if key not in cfgdict:
                    cfgdict[key] = defaultcfgdict[key]

        for key in cfgdict:
            if key == 'load': 
                self.externs, self.externflags = _Parsing.get_seqlist(cfgdict[key])

            elif key == 'ctitle':
                self.ctitle = cfgdict.get('ctitle')
            elif key == 'legend':
                self.legend = cfgdict.getboolean('legend')
            elif key == 'ltitle':
                self.ltitle = cfgdict.get('ltitle')
            elif key == 'rtitle':
                self.rtitle = cfgdict.get('rtitle')
            elif key == 'title':
                self.title = cfgdict.get('title')
            elif key == 'xlabel':
                self.xlabel = cfgdict.get('xlabel')
                self._lock_xlabel = True
            elif key == 'xlim':
                self.xlim = [float(x) for x in _Parsing.get_list(cfgdict.get('xlim', ''))]
                self._lock_xlim = True
            elif key == 'xscale':
                self.xscale = cfgdict.getfloat('xscale')
            elif key == 'ylabel':
                self.ylabel = cfgdict.get('ylabel')
                self._lock_ylabel = True
            elif key == 'ylim':
                self.ylim = [float(x) for x in _Parsing.get_list(cfgdict.get('ylim', ''))]
                self._lock_ylim = True
            elif key == 'yscale':
                self.yscale = cfgdict.getfloat('yscale')
            elif key == 'length':
                self.length = cfgdict.getfloat('length')
                if self.length: self.length = float(self.length)
                self._lock_length = True

            elif key.startswith('mode'):
                self.modes = _Parsing.get_list(cfgdict.get(key))
            elif key.startswith('kernel'):
                self.kernels = _Parsing.get_list(cfgdict.get(key))
            elif key.startswith('addregion'):
                regstr = cfgdict.get(key)
                for obj in _Parsing.get_regions(regstr):
                    #need a better way to hash things
                    #maybe track class and index?
                    for i in range(100, 100+len(self.entities)+1):
                        if 'region{}'.format(i) not in self.entities:
                            self.entities['region{}'.format(i)] = obj
                            break
                #self.modes = _Parsing.get_list(cfgdict.get('mode'))

            elif key.startswith('addtms'):
                tmstr = cfgdict.get(key)
                for obj in _Parsing.get_tmss(tmstr):
                    #need a better way to hash things
                    #maybe track class and index?
                    for i in range(100, 100+len(self.entities)+1):
                        if 'hmmtop{}'.format(i) not in self.entities:
                            self.entities['hmmtop{}'.format(i)] = obj
                            break
            elif key.startswith('notms'):
                if cfgdict.get(key).lower() == 'all': self.notms = SetOfAllThings()
                else: self.notms = _Parsing.get_list(cfgdict.get(key), dtype=int)

            elif key.startswith('addwall'):
                wallstr = cfgdict.get(key)
                for obj in _Parsing.get_walls(wallstr):
                    for i in range(200, 200+len(self.entities)+1):
                        if 'wall{}'.format(i) not in self.entities:
                            self.entities['wall{}'.format(i)] = obj
                            break
            elif '.' in key:
                entname, propname = key.split('.')
                if propname == 'font': self.add_prop(entname, propname, cfgdict.getfloat(key))
                elif propname == 'color': self.add_prop(entname, propname, cfgdict.get(key))
                elif propname == 'linewidth': self.add_prop(entname, propname, cfgdict.getfloat(key))
                elif 'size' in propname: self.add_prop(entname, propname, cfgdict.getfloat(key))
                elif 'alpha' in propname: self.add_prop(entname, propname, cfgdict.getfloat(key))
                else: self.add_prop(entname, propname, cfgdict.get(key))

                #elif not STRICT: warn('Unrecognized key in subplots: "{}"'.format(key))
                #else: error('Unrecognized key in subplots: "{}"'.format(key))

            elif not STRICT: warn('Unrecognized key in subplots: "{}"'.format(key))
            else: error('Unrecognized key in subplots: "{}"'.format(key))


        self.modes = _Parsing.get_list(cfgdict.get('mode', ''))
        self.kernels = _Parsing.get_list(cfgdict.get('kernel', 'flat'))
        self.windows = [int(x) for x in _Parsing.get_list(cfgdict.get('window', '19'))]
        self.multiple = cfgdict.get('multiple', 'overlay')

        #do stuff

    def add_prop(self, entname, propname, value, force=True):
        if force:
            if entname not in self.props: self.props[entname] = {}
            if propname not in self.props[entname]: self.props[entname][propname] = {}
            self.props[entname][propname] = value
        else:
            try: self.props[entname][propname]
            except KeyError: self.add_prop(entname, propname, value)

    def render(self):
        if self.multiple.startswith('frag'): 
            allsequences = []
            for i, (fn, flag) in enumerate(self.externs, self.externflags):
                batch = _Parsing.load_sequences(fn)
                allsequences.extend(batch)
                sequenceflags.extend([flag] * len(batch))

            sequences = []
            fragsequences = []
            for i, seq in enumerate(allsequences):
                if i % 2: sequences.append(seq)
                else: fragsequences.append(seq)

            correction = 0

            if len(self.modes) == 0:
                self.modes = (['hydropathy'] * len(sequences)) + (['hmmtop'] * len(sequences))

            trueseq = len(sequences)
            i = 0
            while len(sequences) < len(self.modes):
                sequences.append(sequences[i])
                fragsequences.append(fragsequences[i])
                i += 1

            for seqi, seq in enumerate(sequences):
                modestr = self.modes[(seqi + correction) % len(self.modes)]
                entname = '{}{}'.format(modestr, seqi)

                mode = libquod.entities.Parser.parse_mode(modestr)
                kernel = self.kernels[(seqi + correction) % len(self.kernels)]
                window = self.windows[(seqi + correction) % len(self.windows)]

                notms = True if (modestr == 'hmmtop') and (seqi in self.notms) else None
                obj = mode.compute(seq, window=window, kernel=kernel, fragment=fragsequences[seqi], notms=notms)
                obj.edgecolor = libquod.entities.get_darkcolor(seqi % trueseq)
                obj.facecolor = libquod.entities.get_lightcolor(seqi % trueseq)

                self.entities[entname] = obj

                
                newxlim, newylim = libquod.entities.update_lims(self.xlim, self.ylim, obj.get_bounding_box())
                if not self._lock_xlim: self.xlim = newxlim
                if not self._lock_ylim: self.ylim = newylim

            
        elif self.multiple.startswith('msa'): pass

        #elif self.multiple == 'overlay':
        #Just assume it's an overlay unless otherwise specified
        else:
            sequences = []
            flags = []
            truesequencecount = 0
            for i, (fn, flag) in enumerate(zip(self.externs, self.externflags)):
                batch = list(_Parsing.load_sequences(fn))
                sequences.extend(batch)
                flags.extend([flag] * len(batch))
                truesequencecount += sum([1 for x in flags if x != 2])

            correction = 0

            if len(self.modes) == 0:
                self.modes = (['hydropathy'] * truesequencecount) + (['hmmtop'] * truesequencecount)
                #self.modes = ['hydropathy', 'hmmtop'] * len(sequences)

            trueseq = len(sequences)
            i = 0
            while len(sequences) < len(self.modes):
                sequences.append(sequences[i])
                flags.append(flags[i])
                i += 1

            for seqi, (seq, flag) in enumerate(zip(sequences, flags)):

                if flag == 2: continue

                modestr = self.modes[(seqi + correction) % len(self.modes)]
                entname = '{}{}'.format(modestr, seqi)
                mode = libquod.entities.Parser.parse_mode(modestr)
                kernel = self.kernels[(seqi + correction) % len(self.kernels)]
                window = self.windows[(seqi + correction) % len(self.windows)]

                notms = True if (modestr == 'hmmtop') and (seqi in self.notms) else None

                if flag == 1:
                    obj = mode.compute(sequences[seqi+1], window=window, kernel=kernel, notms=notms, fragment=seq)
                else:
                    obj = mode.compute(seq, window=window, kernel=kernel, notms=notms)

                obj.edgecolor = libquod.entities.get_darkcolor(seqi % trueseq)
                obj.facecolor = libquod.entities.get_lightcolor(seqi % trueseq)

                self.entities[entname] = obj

                for propent in self.props:
                    if propent == 'all': 
                        for entname in self.entities:
                            for propname in self.props[propent]:
                                if propname == 'color':
                                    setattr(obj, 'edgecolor', self.props[propent][propname])
                                    setattr(obj, 'facecolor', self.props[propent][propname])
                                else:
                                    setattr(obj, propname, self.props[propent][propname])
                    else: 
                        matches = self._fuzzy_find(propent)
                        if matches:
                            for k in matches:
                                for propname in self.props[k]: 
                                    if propname == 'color':
                                        setattr(obj, 'edgecolor', self.props[k][propname])
                                        setattr(obj, 'facecolor', self.props[k][propname])
                                    else:
                                        setattr(obj, propname, self.props[k][propname])
                        else:
                            pass
                            #print(self.entities, matches)

                
                newxlim, newylim = libquod.entities.update_lims(self.xlim, self.ylim, obj.get_bounding_box())
                if not self._lock_xlim: self.xlim = newxlim
                if not self._lock_ylim: self.ylim = newylim


        if not self._lock_length: self.length = self.xlim[1] - self.xlim[0]

        for propent in self.props:
            matches = self._fuzzy_find(propent)
            if matches:
                for k in matches:
                    for propname in self.props[propent]: 
                        setattr(self.entities[k], propname, self.props[propent][propname])
            elif propent in IMPLICIT_ENTITIES: pass
            else: warn('Could not find entity {} (props: {})'.format(propent, self.props[propent]))
    

    def _fuzzy_find(self, q):
        if '-' in q: raise NotImplementedError('Ranged selection not implemented')
        #specific entid
        elif re.search('[0-9]+$', q):
            if q in self.entities: return [q]
            else: 
                out = []
                label = re.findall('^[^-0-9,]+', q)[0]
                rangestr = re.findall('[-0-9,]+$', q)[0]

                #indices = libquod.indexing.unroll_ranges(libquod.entities.Parser.expand_ranges(libquod.entities.Parser.parse_ranges(rangestr, allow_points=True, dtype=int)))
                indices = libquod.entities.Parser.expand_ranges(libquod.entities.Parser.parse_ranges(rangestr, allow_points=True, dtype=int))
                
                for entname in self.entities:
                    if isinstance(entname, str) and entname.startswith(label): 
                        entidlist = re.findall('[0-9]+$', entname)
                        for i in indices:
                            if entidlist and entidlist[0] == str(i): out.append(entname)
                    elif entname == label: out.append(entname)
        #all ents of a type
        else:
            out = []
            label = re.findall('^[^-0-9,]+', q)[0]
            for entname in self.entities:
                if isinstance(entname, str) and entname.startswith(label): out.append(entname)
                elif entname == label: out.append(entname)
        return out
            
            

class AddParser(argparse.ArgumentParser):
    def exit(sel, status=0, message=None): pass #maybe send a warning instead

if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('infile', help='Config file')
    parser.add_argument('--show', action='store_true', help='Display the resulting figure')
    parser.add_argument('--dump', action='store_true', help='Dump entity names')

    args = parser.parse_args()

    cfg = MultiquodConfig.load(args.infile)
    cfg.render()
    cfg.draw(quiet=not args.show)

    if args.dump:
        for subplot in cfg.subplots:
            print('SUBPLOT:', subplot)
            print(','.join([str(entity) for entity in cfg.subplots[subplot].entities]))
