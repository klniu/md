#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
##############################################################################
# @author Xiang Gao
# @email hughgao01@gmail.com
# @Licence GPL v2
##############################################################################
import os
from decimal import Decimal
import argparse
import moltoolkit

def main():
    parser = argparse.ArgumentParser(description='Reassign the atoms charge according the symmetrys of the molecule.\n\nThe format of charge file must be like:\nindex\tcharge\nindex\tcharge\n...', formatter_class=argparse.RawTextHelpFormatter)
    parser.add_argument('-m', '--molfile', required=True, help='Molecule file')
    group = parser.add_mutually_exclusive_group()
    group.add_argument('-s', '--symmetrys', help='Symmetrys, optional. Please give a string whose format is python list. e.g. [[1, 2, 3], [4, 5]]. If you do not assign, the program will calculate it.')
    group.add_argument('--extrasymmetrys', help='Extra Symmetrys, optional. Some symmetrys like the oxygen atoms in sulfonate will not be indentied by program, but you can specify them using this option. Please give a string whose format is python list. e.g. [[1, 2, 3], [4, 5]].')
    parser.add_argument('-c', '--chargefile', required=True, help='Charges file')
    parser.add_argument('-e', '--charge', type=int, default=0, help='Charge of the whole molecule')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    args = parser.parse_args()

    mol = moltoolkit.Mol(args.molfile)
    # Get molecule symmetrys
    if args.symmetries:
        symmetrys = eval(args.symmetry)
    else:
        symmetrys = mol.symmetries

    if args.extrasymmetrys:
        symmetrys += eval(args.extrasymmetrys)
    # symmetrys = [[1, 7], [2, 3, 4, 8, 9, 10]]
    if len(symmetrys) == 0:
        print('Error: Get symmetrys of molecule failure.')
        exit(1)

    # Analyze charges file
    if not os.path.exists(args.chargefile):
        print('Error: Charges file does not exist.')
        exit(1)
    atoms = {}
    # atom index
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
    if len(atoms) != mol.num_atoms:
        print('Error: The atoms number', len(atoms), 'is not equal to the one', mol.num_atoms, 'in the charges file.')
        exit(1)

    # reassign charges
    for sym in symmetrys:
        charge_average = sum([atoms[j] for j in sym]) / len(sym)
        for index in sym:
            try:
                atoms[index] = round(Decimal(charge_average), 3)
            except(Decimal.InvalidOperation):
                print("Fatal Error: The format of the charge file is incorrect.")
                exit(1)

    while True:
        # charge groups
        heavy_charges = {}
        all_atoms_charge = Decimal(0.0)
        for atom in mol.atoms:
            all_atoms_charge += atoms[atom.idx]
            if not atom.is_hydrogen():
                idx = atom.idx
                heavy_charges[idx] = atoms[idx]

        # Split charges into positive and negative
        positive = {}
        negative = {}
        for i in list(heavy_charges.keys()):
            if heavy_charges[i] >= 0:
                positive[i] = heavy_charges[i]
            else:
                negative[i] = heavy_charges[i]

        # record the max or min charge
        max_key = max(positive.values())
        max_keys = [i for i in positive.keys() if positive[i] == max_key]
        min_key = min(negative.values())
        min_keys = [i for i in negative.keys() if negative[i] == min_key]

        spare_charge = all_atoms_charge - Decimal(args.charge)
        if spare_charge == 0:
            break
        if spare_charge > 0:
            for i in max_keys:
                atoms[i] -= round(spare_charge / len(max_keys), 3)
        else:
            for i in min_keys:
                atoms[i] -= round(spare_charge / len(min_keys), 3)

    # Check if the charge of the whole molecule is integer
    charges = sum(atoms)
    if round(charges, 3) != round(charges, 0):
        print("Error: The charges of the whole molecule is not a interger. Something is wrong.")
        exit(1)

    with open(args.output, 'w') as f:
        for i in sorted(atoms.keys()):
            f.write(str(atoms[i]).rjust(6) + '\n')

if __name__ == '__main__':
    main()
