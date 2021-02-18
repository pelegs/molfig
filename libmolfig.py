#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

from colorsys import rgb_to_hls, hls_to_rgb
#import drawSvg as draw
import svgwrite
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

def format_color(color):
    return 'rgb({})'.format(','.join(map(str, color)))

def create_atom(drawing, element='C', pos=(0,0)):
    px, py = pos
    radius = atoms[element]['radius']
    ed = atoms[element]
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

    # effects
    blur = drawing.defs.add(drawing.filter())
    blur.feGaussianBlur(in_='SourceGraphic', stdDeviation=3)
    blur_filter = blur.get_funciri()

    # shapes
    main_circle = drawing.circle(pos, radius, fill=main_gradient_paintsever)
    drawing.add(main_circle)
    center_blur = drawing.ellipse(pos, (radius*0.6, radius*0.4),
                                  fill=brighter_color, fill_opacity=0.6,
                                  filter=blur_filter)
    drawing.add(center_blur)
    upper_glow = drawing.ellipse((px, py-radius*0.37), (radius*0.8, radius*0.55),
                                 fill=upper_glow_gradient_paintsever)
    drawing.add(upper_glow)

    # label
    #label = draw.Text(element, x=px, y=py, text_anchor='middle',
    #                  font_family='FreeSans', font_weight='bold', fontSize=20)
    #drawing.append(label)


if __name__ == '__main__':
    d = svgwrite.Drawing(filename='images/first_test.svg')
    for key in atoms:
        create_atom(d, key, np.random.uniform((-250,-250), (250,250), 2))
    d.save()
