#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

from colorsys import rgb_to_hls, hls_to_rgb
#import drawSvg as draw
import svgwrite
import numpy as np
from subprocess import call
import openbabel
from openbabel import pybel


atoms = {
        1:  {'symbol': 'H',  'radius': 15, 'col': (170, 204, 255), 'text_col': 'black'},
        4:  {'symbol': 'B',  'radius': 17, 'col': (123, 76, 0),    'text_col': 'white'},
        6:  {'symbol': 'C',  'radius': 19, 'col': (26, 26, 26),    'text_col': 'white'},
        7:  {'symbol': 'N',  'radius': 19, 'col': (0, 73, 228),    'text_col': 'white'},
        8:  {'symbol': 'O',  'radius': 20, 'col': (212, 0, 0),     'text_col': 'white'},
        9:  {'symbol': 'F',  'radius': 18, 'col': (34, 193, 71),   'text_col': 'white'},
        14: {'symbol': 'Si', 'radius': 20, 'col': (86, 83, 86),    'text_col': 'white'},
        15: {'symbol': 'P',  'radius': 19, 'col': (36, 207, 128),  'text_col': 'black'},
        16: {'symbol': 'S',  'radius': 20, 'col': (240, 209, 0),   'text_col': 'black'},
        17: {'symbol': 'Cl', 'radius': 30, 'col': (163, 227, 0),   'text_col': 'black'},
        33: {'symbol': 'As', 'radius': 19, 'col': (199, 90, 152),  'text_col': 'black'},
        34: {'symbol': 'Se', 'radius': 19, 'col': (97, 208, 255),  'text_col': 'black'},
        35: {'symbol': 'Br', 'radius': 31, 'col': (84, 6, 6),      'text_col': 'white'},
        53: {'symbol': 'I',  'radius': 33, 'col': (109, 32, 136),  'text_col': 'white'},
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

def format_color(color):
    return 'rgb({})'.format(','.join(map(str, color)))

def create_lower_glow(center, radius):
    shape_start = 'm {},{} '.format(center[0]+radius*0.7, center[1])
    shape_coords = np.array([
            [0.0,-0.382667,-0.329,0.270667,-0.714,0.270667],
            [-0.385,0.0,-0.681333,-0.653333,-0.681333,-0.270667],
            [0.0,0.382667,0.312667,0.693,0.697667,0.693],
            [0.385,0.0,0.697667,-0.310333,0.697667,-0.693],
        ]) * radius
    rows = 'c '
    rows += ' '.join(['{},{} {},{} {},{}'.format(row[0], row[1], row[2], row[3], row[4], row[5])
             for row in shape_coords])
    rows += ' z'
    return shape_start + rows

def draw_atom(drawing, atomic_num=6, pos=(0,0), scale=1):
    scaled_pos = np.array(pos) * scale
    px, py = scaled_pos
    symbol = atoms[atomic_num]['symbol']
    radius = atoms[atomic_num]['radius']
    ed = atoms[atomic_num]
    main_color = format_color(ed['col'])
    brighter_color = format_color(lighten_color(ed['col'], 0.8))

    # gradients
    main_gradient = drawing.linearGradient((0,0), (0,1))
    main_gradient.add_stop_color(0, brighter_color)
    main_gradient.add_stop_color(1, main_color)
    drawing.defs.add(main_gradient)
    main_gradient_paintsever = main_gradient.get_paint_server(default='currentColor')

    upper_glow_gradient = drawing.linearGradient((0,0), (0,1))
    upper_glow_gradient.add_stop_color(0, 'white')
    upper_glow_gradient.add_stop_color(0.75, 'white', 0)
    drawing.defs.add(upper_glow_gradient)
    upper_glow_gradient_paintsever = upper_glow_gradient.get_paint_server(default='currentColor')
    
    lower_glow_gradient = drawing.radialGradient((0.5,0.5), 1)
    lower_glow_gradient.add_stop_color(0, 'white')
    lower_glow_gradient.add_stop_color(0.75, 'white', 0)
    drawing.defs.add(lower_glow_gradient)
    lower_glow_gradient_paintsever = lower_glow_gradient.get_paint_server(default='currentColor')

    # effects
    center_blur = drawing.defs.add(drawing.filter())
    center_blur.feGaussianBlur(in_='SourceGraphic', stdDeviation=3)
    center_blur_filter = center_blur.get_funciri()
    
    lower_blur = drawing.defs.add(drawing.filter())
    lower_blur.feGaussianBlur(in_='SourceGraphic', stdDeviation=8)
    lower_blur_filter = lower_blur.get_funciri()

    # shapes
    main_circle = drawing.circle(scaled_pos, radius, fill=main_gradient_paintsever)
    drawing.add(main_circle)
    center_ellipse = drawing.ellipse(scaled_pos, (radius*0.6, radius*0.4),
                                  fill=brighter_color, fill_opacity=0.6,
                                  filter=center_blur_filter)
    drawing.add(center_ellipse)

    upper_glow = drawing.ellipse((px, py-radius*0.37), (radius*0.8, radius*0.55),
                                 fill=upper_glow_gradient_paintsever)
    drawing.add(upper_glow)

    d = create_lower_glow(scaled_pos, radius)
    lower_glow = drawing.path(d=d, fill=lower_glow_gradient_paintsever,
                              filter=lower_blur_filter)
    drawing.add(lower_glow)

    label = drawing.text(
        symbol,
        scaled_pos,
        style='text-anchor:middle;\
               dominant-baseline:middle;\
               font-family:FreeSans;\
               font-weight:bold',
        font_size='22px',
        fill=ed['text_col'],
    )
    drawing.add(label)


if __name__ == '__main__':
    d = svgwrite.Drawing(filename='images/first_test.svg')
    #create_atom(d, 1)
    #for key in atoms:
    #    create_atom(d, key, np.random.uniform((-250,-250), (250,250), 2))

    mol = pybel.readstring('smi', 'CC(=O)C')
    mol.addh()
    mol.make2D()
    for atom in mol.atoms:
        atomic_num = atom.atomicnum
        pos = atom.coords[:2]
        draw_atom(d, atomic_num, pos, scale=75)
    d.save()
