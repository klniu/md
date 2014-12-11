#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
##############################################################################
# @file unitedcharge.py
#
# @date 2012-04-20
# @author Xiang Gao
# @email email@Klniu.com
#
# @Licence GPL v2
#
# @brief
#
# @detail
#
##############################################################################
import os
from decimal import Decimal
import argparse
import pybel


def main():
    parser = argparse.ArgumentParser(description='Convert all atom charges to united charges and renumber the atom index.\n\nThe format of charge file must be like:\nindex\tcharge\nindex\tcharge\n...', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-f', '--format', default='pdb', help='The format of molecule(default=pdb)')
    parser.add_argument('-m', '--molfile', required=True, help='Molecule file')
    parser.add_argument('-c', '--chargefile', required=True, help='Charges file')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    args = parser.parse_args()

    mol = next(pybel.readfile(args.format, args.molfile))
    ob = pybel.ob
    obmol = mol.OBMol
    atomsNum = len(mol.atoms)

    # Analyze charges file
    if not os.path.exists(args.chargefile):
        print('Error: Charges file does not exist.')
        exit(1)
    atoms = {}
    # 原子索引
    index = 0
    with open(args.chargefile) as f:
        for line in f:
            charge = line.strip('\n')
            if charge[0] == '#':
                continue
            index += 1
            try:
                atoms[index] = round(Decimal(charge), 3)
            except(Decimal.InvalidOperation):
                print("Fatal Error: The format of the charge file is incorrect.")
                exit(1)
    if len(atoms) != atomsNum:
        print('Error: The atoms number', len(atoms), 'is not equal to the one', atomsNum, 'in the charges file.')
        exit(1)

    # 检查电荷是否为整数
    charges = sum(atoms)
    if round(charges, 3) != charges:
        print("Error: The charges of the whole molecule is not a interger.")
        exit(1)

    for obatom in ob.OBMolAtomIter(obmol):
        if obatom.GetHyb() == 3 and obatom.GetAtomicNum() == 6:
            for nbr in ob.OBAtomAtomIter(obatom):
                if nbr.IsHydrogen():
                    atoms[obatom.GetIdx()] += atoms[nbr.GetIdx()]
                    atoms.pop(nbr.GetIdx())

    new_idx = 1
    united_mol = {}
    for i in sorted(atoms.keys()):
        united_mol[new_idx] = atoms[i]
        new_idx += 1

    with open(args.output, 'w') as f:
        for i in sorted(united_mol.keys()):
            f.write(str(united_mol[i]).rjust(6) + '\n')

if __name__ == '__main__':
    main()
