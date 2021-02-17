#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

from colorsys import rgb_to_hls, hls_to_rgb
import drawSvg as draw
import numpy as np
from subprocess import call


atoms = {
    'H':  {'radius': 15, 'col': (170, 204, 255), 'text_col': 'black'},
    'B':  {'radius': 17, 'col': (123, 76, 0),    'text_col': 'white'},
    'C':  {'radius': 19, 'col': (26, 26, 26),    'text_col': 'white'},
    'N':  {'radius': 19, 'col': (0, 73, 228),    'text_col': 'white'},
    'O':  {'radius': 20, 'col': (212, 0, 0),     'text_col': 'white'},
    'F':  {'radius': 18, 'col': (34, 193, 71),   'text_col': 'white'},
    'Si': {'radius': 20, 'col': (86, 83, 86),    'text_col': 'white'},
    'P':  {'radius': 19, 'col': (36, 207, 128),  'text_col': 'black'},
    'S':  {'radius': 20, 'col': (240, 209, 0),   'text_col': 'black'},
    'Cl': {'radius': 30, 'col': (163, 227, 0),   'text_col': 'black'},
    'As': {'radius': 19, 'col': (199, 90, 152),  'text_col': 'black'},
    'Se': {'radius': 19, 'col': (97, 208, 255),  'text_col': 'black'},
    'Br': {'radius': 31, 'col': (84, 6, 6),      'text_col': 'white'},
    'I':  {'radius': 33, 'col': (109, 32, 136),  'text_col': 'white'},
}

def rgbhex(c):
    return '#%02x%02x%02x' % c

def adjust_color_lightness(color, factor):
    r, g, b = color
    h, l, s = rgb_to_hls(r / 255.0, g / 255.0, b / 255.0)
    l = max(min(l * factor, 1.0), 0.0)
    r, g, b = hls_to_rgb(h, l, s)
    return int(r * 255), int(g * 255), int(b * 255)

def lighten_color(color, factor=0.1):
    return adjust_color_lightness(color, 1+factor)

def darken_color(color, factor=0.1):
    return adjust_color_lightness(color, 1-factor)

def create_atom(drawing, element='C', pos=(0,0)):
    px, py = pos
    rad = atoms[element]['radius']
    ed = atoms[element]
    main_color = rgbhex(ed['col'])
    brighter_color = rgbhex(lighten_color(ed['col'], 0.8))

    # gradients
    main_gradient = draw.LinearGradient(px+0, py-rad/1.5, px+0, py+rad/1.5)
    main_gradient.addStop(0, main_color, 1)
    main_gradient.addStop(1, brighter_color, 1)

    # shapes
    main_circle = draw.Circle(px, py, atoms[element]['radius'], fill=main_gradient)
    drawing.append(main_circle)

    # label
    label = draw.Text(element, x=px, y=py, text_anchor='middle',
                      font_family='FreeSans', font_weight='bold', fontSize=20)
    drawing.append(label)


if __name__ == '__main__':
    d = draw.Drawing(500, 500, origin='center', displayInline=False)
    for key in atoms:
        create_atom(d, key, np.random.uniform((-250,-250), (250,250), 2))
    d.saveSvg('images/first_test.svg')
