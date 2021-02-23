#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

from colorsys import rgb_to_hls, hls_to_rgb
#import drawSvg as draw
import svgwrite
import numpy as np
from subprocess import call
import openbabel
from openbabel import pybel
from sys import argv
import os


def rotMat(angle):
    c = np.cos(angle)
    s = np.sin(angle)
    return np.array([[c,-s],
                     [s,c]])


atomic_data = {
        1:  {'symbol': 'H',  'radius': 15,
             'color': (170, 204, 255), 'text_color': 'black'},
        4:  {'symbol': 'B',  'radius': 17,
             'color': (123, 76, 0),    'text_color': 'white'},
        6:  {'symbol': 'C',  'radius': 19,
             'color': (26, 26, 26),    'text_color': 'white'},
        7:  {'symbol': 'N',  'radius': 19,
             'color': (0, 73, 228),    'text_color': 'white'},
        8:  {'symbol': 'O',  'radius': 20,
             'color': (212, 0, 0),     'text_color': 'white'},
        9:  {'symbol': 'F',  'radius': 18,
             'color': (34, 193, 71),   'text_color': 'white'},
        14: {'symbol': 'Si', 'radius': 20,
             'color': (86, 83, 86),    'text_color': 'white'},
        15: {'symbol': 'P',  'radius': 19,
             'color': (36, 207, 128),  'text_color': 'black'},
        16: {'symbol': 'S',  'radius': 20,
             'color': (240, 209, 0),   'text_color': 'black'},
        17: {'symbol': 'Cl', 'radius': 30,
             'color': (163, 227, 0),   'text_color': 'black'},
        33: {'symbol': 'As', 'radius': 19,
             'color': (199, 90, 152),  'text_color': 'black'},
        34: {'symbol': 'Se', 'radius': 19,
             'color': (97, 208, 255),  'text_color': 'black'},
        35: {'symbol': 'Br', 'radius': 31,
             'color': (84, 6, 6),      'text_color': 'white'},
        53: {'symbol': 'I',  'radius': 33,
             'color': (109, 32, 136),  'text_color': 'white'},
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


class Atom:
    def __init__(self, atomic_num=6, pos=(0,0), scale=1):
        self.atomic_num = atomic_num
        self.symbol = atomic_data[atomic_num]['symbol']
        self.radius = atomic_data[atomic_num]['radius']
        self.color = atomic_data[atomic_num]['color']
        self.text_color = atomic_data[atomic_num]['text_color']
        self.pos = np.array(pos) * scale

    def create_lower_glow(self):
        shape_start = 'm {},{} '.format(self.pos[0]+self.radius*0.7,
                                        self.pos[1])
        shape_coords = np.array([
            [0.0,-0.382667,-0.329,0.270667,-0.714,0.270667],
            [-0.385,0.0,-0.681333,-0.653333,-0.681333,-0.270667],
            [0.0,0.382667,0.312667,0.693,0.697667,0.693],
            [0.385,0.0,0.697667,-0.310333,0.697667,-0.693],
            ]) * self.radius
        rows = 'c '
        rows += ' '.join(['{},{} {},{} {},{}'.format(row[0], row[1],
                                                     row[2], row[3],
                                                     row[4], row[5])
            for row in shape_coords])
        rows += ' z'
        return shape_start + rows

    def draw(self, drawing):
        px, py = self.pos
        main_color = format_color(self.color)
        brighter_color = format_color(lighten_color(self.color, 0.8))

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
        main_circle = drawing.circle(self.pos, self.radius, fill=main_gradient_paintsever)
        drawing.add(main_circle)
        center_ellipse = drawing.ellipse(self.pos, (self.radius*0.6, self.radius*0.4),
                fill=brighter_color, fill_opacity=0.6,
                filter=center_blur_filter)
        drawing.add(center_ellipse)

        upper_glow = drawing.ellipse((px, py-self.radius*0.37),
                                     (self.radius*0.8, self.radius*0.55),
                fill=upper_glow_gradient_paintsever)
        drawing.add(upper_glow)

        d = self.create_lower_glow()
        lower_glow = drawing.path(d=d, fill=lower_glow_gradient_paintsever,
                filter=lower_blur_filter)
        drawing.add(lower_glow)

        label = drawing.text(
                self.symbol,
                self.pos,
                style='text-anchor:middle;\
                        dominant-baseline:central;\
                        font-family:FreeSans;\
                        font-weight:bold',
                        font_size='22px',
                        fill=self.text_color,
                        )
        drawing.add(label)

    def __str__(self):
        return '{}, ({},{})'.format(self.atomic_num, self.pos[0], self.pos[1])


class Bond:
    def __init__(self, type=1, startAtom=None, endAtom=None):
        self.type = type
        self.startAtom = startAtom
        self.endAtom = endAtom

    def __str__(self):
        #return '{} {} {}'.format(self.type, self.startIdx, self.endIdx)
        return str(self.type)

    def draw(self, drawing):
        if self.type == 1:
            line = drawing.line(start=self.startAtom.pos,
                                end=self.endAtom.pos,
                                stroke='black',
                                stroke_width=7)
            drawing.add(line)
        elif self.type == 2:
            start = self.startAtom.pos
            end = self.endAtom.pos
            direct = end - start
            direct_unit = direct / np.linalg.norm(direct)
            perp = np.array([direct_unit[1], -direct_unit[0]])
            s1 = start + 7*perp
            s2 = start - 7*perp
            e1 = end + 7*perp
            e2 = end - 7*perp
            line1 = drawing.line(start=s1,
                                 end=e1,
                                 stroke='black',
                                 stroke_width=7)
            line2 = drawing.line(start=s2,
                                 end=e2,
                                 stroke='black',
                                 stroke_width=7)
            drawing.add(line1)
            drawing.add(line2)
        else:
            pass


class Molecule:
    def __init__(self, atoms=[], bonds=[], scale=1):
        self.atoms = atoms
        self.bonds = bonds
        self.scale = scale

    def from_OBMol(self, mol):
        # Create atoms
        for atom in mol.atoms:
            self.atoms.append(Atom(atomic_num=atom.atomicnum,
                                   pos=atom.coords[:2],
                                   scale=self.scale))

        # Create bonds
        mol.OBMol.AddHydrogens()
        self.num_bonds = mol.OBMol.NumBonds()
        for i in range(self.num_bonds):
            bond = mol.OBMol.GetBond(i)
            si = bond.GetBeginAtomIdx()-1
            ei = bond.GetEndAtomIdx()-1
            self.bonds.append(Bond(bond.GetBondOrder(),
                                   self.atoms[si],
                                   self.atoms[ei]))

    def __str__(self):
        s  = '\n'.join([atom.__str__() for atom in self.atoms])
        s += '\n'.join([bond.__str__() for bond in self.bonds])
        return s

    def draw(self, drawing):
        for bond in self.bonds:
            bond.draw(drawing)
        for atom in self.atoms:
            atom.draw(drawing)

    def rotate(self, angle=0):
        mat = rotMat(angle)
        for atom in self.atoms:
            relPos = atom.pos - self.center
            new_relPos = np.dot(mat, relPos)
            atom.pos = new_relPos + self.center

    def create_bounding_box(self):
        self.xmin = min([atom.pos[0]-atom.radius for atom in self.atoms])
        self.xmax = max([atom.pos[0]+atom.radius for atom in self.atoms])
        self.ymin = min([atom.pos[1]-atom.radius for atom in self.atoms])
        self.ymax = max([atom.pos[1]+atom.radius for atom in self.atoms])
        self.bounding_box = ((self.xmin, self.ymin), (self.xmax, self.ymax))
        self.center = np.array([self.xmin+(self.xmax-self.xmin)/2,
                                self.ymin+(self.ymax-self.ymin)/2])

    def position_center(self, new_center):
        trans = new_center - self.center
        for atom in self.atoms:
            atom.pos += trans
        self.create_bounding_box()


if __name__ == '__main__':
    SMILES = argv[1]
    SCALE = int(argv[2])

    mol = pybel.readstring('smi', SMILES)
    mol.addh()
    mol.make2D()

    molFigure = Molecule(scale=SCALE)
    molFigure.from_OBMol(mol)
    molFigure.create_bounding_box()
    #width, height = (molFigure.xmax-molFigure.xmin, molFigure.ymax-molFigure.ymin)
    width, height = 1200, 1200
    molFigure.position_center(np.array([width/2, height/2]))

    d = svgwrite.Drawing(filename='video/frame{:03d}.svg'.format(frame),
                         size=('{}px'.format(width), '{}px'.format(height)))
    molFigure.draw(d)
    d.save()
