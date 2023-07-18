from . import entities
from . import indexing

import re
import subprocess
import numpy as np
import matplotlib.cm
import sys

def info(*things):
    print('[INFO]', *things, file=sys.stderr)
def warn(*things):
    print('[WARNING]', *things, file=sys.stderr)
def error(*things):
    print('[ERROR]', *things, file=sys.stderr)
    exit(1)

defaults = {
    "yspan":[-2.75, -2.5],
}
def draw_other(ax, **kwargs):
    xlim = kwargs.get('xlim')
    ylim = kwargs.get('ylim')
    if 'addtms' in kwargs and kwargs['addtms']:
        for batch in kwargs['addtms']:
            args = batch.split(':')
            if len(args) == 0: raise ValueError('Invalid TMS specification: "{}"'.format(batch))

            if len(args) >= 1: spans = entities.Parser.parse_ranges(args[0])
            else: spans = []

            if len(args) >= 2: color = args[1]
            else: color = 'orange'

            if len(args) >= 3: alpha = float(args[2])
            else: alpha = None

            if alpha is None: hmmtop = entities.HMMTOP(spans=spans, fc=color)
            else: hmmtop = entities.HMMTOP(spans=spans, fc=color, alpha=alpha)

            hmmtop.plot(ax)
            xlim, ylim = update_lims(xlim, ylim, hmmtop.get_bounding_box())

    if 'addregion' in kwargs and kwargs['addregion']:
        #spans:labeltext:y:facecolor(,textcolor):valign,halign

        for batch in kwargs['addregion']:
            args = batch.split(':')
            spans = entities.Parser.parse_ranges(args.pop(0))

            #incerti ordinis {
            try: text = args.pop(0)
            except IndexError: text = ''

            regionkwargs = {'alpha':1.0}

            try: 
                rawy = args.pop(0)
                y = [float(_) for _ in rawy.split(',')]
                if len(y) == 2:
                    if y[0] == y[-1]:
                        yspan = [y[0]-0.15, y[0]]
                    else:
                        yspan = y
                elif len(y) == 1:
                    yspan = [y[0]-0.15, y[0]]
            except IndexError:
                yspan = [-2.8, -2.8+0.15]

            try:
                colors = args.pop(0)
                if len(colors.split(',')) == 1:
                    regionkwargs['facecolor'] = colors
                else:
                    regionkwargs['facecolor'], regionkwargs['textcolor'] = colors.split(',')
            except IndexError:
                regionkwargs['facecolor'] = 'blue'
                regionkwargs['textcolor'] = 'white'

            #}

            try:
                align = args.pop(0)
                if len(align) == 1:
                    regionkwargs['valign'] = align
                    regionkwargs['halign'] = 'c'
                else: regionkwargs['valign'], regionkwargs['halign'] = align
            except IndexError:
                regionkwargs['valign'] = 'm'
                regionkwargs['halign'] = 'c'


            region = entities.Region(spans=spans, yspan=yspan, text=text, **regionkwargs)
            #region = entities.Region(spans=spans, yspan=yspan, fc=color, text=text, ha='c', va='t', fontsize=kwargs.get('fontsize'))
            region.plot(ax)
            xlim, ylim = update_lims(xlim, ylim,region.get_bounding_box())

    if 'adddashedregion' in kwargs and kwargs['adddashedregion']:
        #spans:labeltext:y:facecolor(,textcolor):valign,halign

        for batch in kwargs['adddashedregion']:
            args = batch.split(':')
            spans = entities.Parser.parse_ranges(args.pop(0))

            #incerti ordinis {
            try: text = args.pop(0)
            except IndexError: text = ''

            regionkwargs = {'alpha':1.0}

            try: 
                rawy = args.pop(0)
                y = [float(_) for _ in rawy.split(',')]
                if len(y) == 2:
                    if y[0] == y[-1]:
                        yspan = [y[0]-0.15, y[0]]
                    else:
                        yspan = y
                elif len(y) == 1:
                    yspan = [y[0]-0.15, y[0]]
            except IndexError:
                yspan = [-2.8, -2.8+0.15]

            try:
                colors = args.pop(0)
                if len(colors.split(',')) == 1:
                    regionkwargs['facecolor'] = colors
                else:
                    regionkwargs['facecolor'], regionkwargs['textcolor'] = colors.split(',')
            except IndexError:
                regionkwargs['facecolor'] = 'blue'
                regionkwargs['textcolor'] = 'white'

            #}

            try:
                align = args.pop(0)
                if len(align) == 1:
                    regionkwargs['valign'] = align
                    regionkwargs['halign'] = 'c'
                else: regionkwargs['valign'], regionkwargs['halign'] = align
            except IndexError:
                regionkwargs['valign'] = 'm'
                regionkwargs['halign'] = 'c'

            regionkwargs["dashed"] = kwargs.get("dashed_fraction")
            regionkwargs["period"] = kwargs.get("dashed_period")


            region = entities.Region(spans=spans, yspan=yspan, text=text, **regionkwargs)
            #region = entities.Region(spans=spans, yspan=yspan, fc=color, text=text, ha='c', va='t', fontsize=kwargs.get('fontsize'))
            region.plot(ax)
            xlim, ylim = update_lims(xlim, ylim,region.get_bounding_box())
    if kwargs.get('addwalls'):
        for batch in kwargs['addwalls']:
            args = batch.split(':')
            spans = entities.Parser.parse_ranges(args[0], allow_points=True)

            #if len(args) >= 2: scale = float(args[1])
            #else: scale = 1.0

            if len(args) >= 2: y = float(args[1])
            else: y = 1.732

            if len(args) >= 3: ypos = args[2]
            else: ypos = '+-'

            if len(args) >= 4: text = args[3]
            else: text = ''

            walls = entities.Wall(spans=spans, y=y, ypos=ypos, text=text)
            xlim, ylim = update_lims(xlim, ylim, walls.get_bounding_box())
            walls.plot(ax)

    if kwargs.get('addwedges'):
        for batch in kwargs['addwedges']:
            args = batch.split(':')
            x = [(float(c), float(c)) for c in args[0].split(",")]

            if len(args) >= 2: y = float(args[1])
            else: y = 1.732

            if len(args) >= 3: ypos = args[2]
            else: ypos = "+-"

            if len(args) >= 4: 
                try: scale = float(args[3])
                except ValueError as e:
                    if args[3] == "+": scale = 1.0
                    elif args[3] == "-": scale = -1.0
                    else: raise e
            else: scale = 1.0

            walls = entities.Wall(spans=x, y=y, ypos=ypos, scale=scale)
            xlim, ylim = update_lims(xlim, ylim, walls.get_bounding_box())
            walls.plot(ax)

            

        
    return xlim, ylim

def onlycaps(text):
    if re.match('^[-A-Z\s]+$', text):
        return True

def update_lims(oldxlim, oldylim, newlims=None):
    #none, none -> none:
    if newlims is None:
        xlim = oldxlim
        ylim = oldylim
    #none, none -> none, none
    elif newlims[0] is None and newlims[1] is None:
        xlim = oldxlim
        ylim = oldylim
    elif newlims[0] is not None or newlims[1] is not None:
        if newlims[0] is None: xlim = oldxlim
        else:
            if oldxlim is None: xlim = newlims[0]
            else: xlim = [min(oldxlim[0], newlims[0][0]), max(oldxlim[1], newlims[0][1])]
        if newlims[1] is None: ylim = oldylim
        else:
            if oldylim is None: ylim = newlims[1]
            else: ylim = [min(oldylim[0], newlims[1][0]), max(oldylim[1], newlims[1][1])]
    else: raise ZeroDivisionError("impossible exception")

    return xlim, ylim

def validate_mode(mode):
    if mode.startswith('hydro'): return True
    elif mode.startswith('amphi'): return True
    elif mode.startswith('charge'): return True
    elif mode.startswith('entropy'): return True
    elif mode.startswith('psipred'): return True
    elif mode.startswith('conform'): return True
    elif mode.startswith('ident'): return True

    else: error('Mode "{}" not implemented'.format(mode))
