#!/usr/bin/env python3
#
# Copyright 2013 Xiang Gao
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#   http://www.apache.org/licenses/LICENSE-2.0
#
'''Remove hydrogens from the molecule.'''

import pybel
import openbabel
import argparse


def main():
    parser = argparse.ArgumentParser(description='Remove hydrogens from the molecule.')
    parser.add_argument('-m', '--mol', action='store', help='Molecule file')
    parser.add_argument('-s', '--save', action='store', help='The file saving the result to')
    args = parser.parse_args()
    mol = next(pybel.readfile(args.mol[-3:], args.mol))
    obmol = mol.OBMol

    test1 = 0
    test2 = 0
    test3 = 0
    finish = False
    deleted = False
    while not finish:
        for obatom in openbabel.OBMolAtomIter(obmol):
            if obatom.IsHydrogen() or not obatom.IsCarbon():
                continue
            if obatom.GetHyb() == 3:
                for neibor_atom in openbabel.OBAtomAtomIter(obatom):
                    if neibor_atom.IsHydrogen():
                        obmol.DeleteAtom(neibor_atom)
                        deleted = True
                        break
                    else:
                        deleted = False
            if deleted:
                break
        else:
            finish = True

    mol.write(args.save[-3:], args.save)

if __name__ == '__main__':
    main()
