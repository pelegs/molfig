#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import openbabel
from openbabel import pybel

mol = pybel.readstring('smi', 'CCO')
mol.addh()
mol.make2D()
for atom in mol.atoms:
    print(atom.atomicnum, atom.coords)
