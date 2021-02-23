#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

from openbabel import openbabel
from openbabel import pybel


SMILES = 'CC(=O)C'


mol = pybel.readstring('smi', SMILES)
mol.addh()
mol.make2D()
for i, atom in enumerate(mol.atoms):
    print(i+1, atom.atomicnum, atom.coords)

mol.OBMol.AddHydrogens()
num_bonds = mol.OBMol.NumBonds()
for i in range(num_bonds):
    bond = mol.OBMol.GetBond(i)
    print(i+1, bond.GetBondOrder(), bond.GetBeginAtomIdx(), bond.GetEndAtomIdx())
