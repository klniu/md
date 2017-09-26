#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# Copyright (C) Hugh Gao
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
##
# @file fragitpjoin.py
# @brief According the indics match of a molecule and its fragment, convert the itp file of fragment to molecule to use.
# @author Hugh Gao, hugh8505@gmail.com
# @version 0.1alpa
# @date 2012-04-21
import moltoolkit
from gromacs.fileformats.itp import ITP
from os.path import basename


def main():
    import argparse
    parser = argparse.ArgumentParser(description='According the indics match of a molecule and its fragment, convert the itp file of fragment to molecule to use. The second molecule should be a fragment or identical with the first one.')
    parser.add_argument('-1', '--mol1', required=True, help='The first Molecule file')
    parser.add_argument('-2', '--mol2', required=True, help='The second molecule file (fragment).')
    parser.add_argument('-t', '--top2', required=True, help='The topology file of the second molecule (fragment).')
    parser.add_argument('-p', '--pairs', help='Optional. The match pairs of atoms indices in the second molecule and the first molecule, optional. Please give a string whose format is python list. e.g. [(1, 2), (4, 5)]. If you do not assign, the program will calculate it. ')
    parser.add_argument('-o', '--output', required=True, help='Output file')
    #parser.add_argument('--include_atoms', help='--include_atoms "[atom_index, ...]", 指定第一个分子中必须包含的原子索引(从1开始)，用于第一个分子中包含多个匹配片段时指定匹配的部分, 例如[1,2]')
    args = parser.parse_args()

    itpfile = args.top2
    # Get matches
    mol = moltoolkit.Mol(args.mol1)
    if args.pairs:
        matches = eval(args.pairs)
    else:
        #if args.include_atoms:
        #    include_atoms = eval(args.include_atoms)
        #    if isinstance(include_atoms, list):
        #        matches = mol.getSubstructureMaps(moltoolkit.Mol(args.mol2), include_atoms)
        #    else:
        #        print('Error: 错误的include_atoms格式，应类似为"[1,2,3]"')
        #        exit(1)
        #else:
        matches = mol.getSubstructureMaps(moltoolkit.Mol(args.mol2))

    if len(matches) == 0:
        print("Error: The small molecule is not a fragment of the molecule.")
        exit(1)
    matches_dict = {}
    for i, j in matches:
        matches_dict[i] = j

    # Parser itp
    itp = ITP(itpfile)

    # [atoms]
    atoms_data = itp.header.moleculetype.atoms.data
    for atom in atoms_data:
        idx = matches_dict[atom[0]]
        # exchange the index
        atom[0] = idx
        # get the atom type to replace the origin ones in itp
        atom[4] = mol.get_atom(idx).atomID
    atoms_data.sort()

    # make cgnr
    for idx, atom in enumerate(atoms_data):
        equal_idx = []
        for next_idx in range(idx, len(atoms_data)):
            if atoms_data[next_idx][5] == atom[5]:
                equal_idx.append(next_idx)
            else:
                break
        if idx == 0 and atom[5] != 1:
            for idx_new in equal_idx:
                atoms_data[idx_new][5] = 1
        if idx != 0 and atom[5] != atoms_data[idx - 1][5]:
            last_cgr = atoms_data[idx - 1][5]
            for idx_new in equal_idx:
                atoms_data[idx_new][5] = last_cgr + 1

    # [bonds]
    for bond in itp.header.moleculetype.bonds.data:
        old_idx0, old_idx1 = bond[0], bond[1]
        idx0, idx1 = matches_dict[old_idx0], matches_dict[old_idx1]
        bond[0] = min(idx0, idx1)
        bond[1] = max(idx0, idx1)
        bond[5] = basename(args.output) + "  " + str(old_idx0) + ' ' + str(old_idx1)
    itp.header.moleculetype.bonds.data.sort()

    # [pairs]
    for pair in itp.header.moleculetype.pairs.data:
        idx0, idx1 = matches_dict[pair[0]], matches_dict[pair[1]]
        pair[0] = min(idx0, idx1)
        pair[1] = max(idx0, idx1)
    itp.header.moleculetype.pairs.data.sort()

    # [angles]
    for angle in itp.header.moleculetype.angles.data:
        old_idx0, old_idx1, old_idx2 = angle[0], angle[1], angle[2]
        idx0, idx1, idx2 = matches_dict[angle[0]], matches_dict[angle[1]], matches_dict[angle[2]]
        angle[0] = min(idx0, idx2)
        angle[1] = idx1
        angle[2] = max(idx0, idx2)
        angle[6] = basename(args.output) + "  " + str(old_idx0) + ' ' + str(old_idx1) + ' ' + str(old_idx2)
    itp.header.moleculetype.angles.set_data(sorted(itp.header.moleculetype.angles.data, key=lambda data: data[1]))

    # [dihedrals]
    for dihedral in itp.header.moleculetype.dihedrals.data:
        old_idx0, old_idx1, old_idx2, old_idx3 = dihedral[0], dihedral[1], dihedral[2], dihedral[3]
        idx0, idx1, idx2, idx3 = matches_dict[dihedral[0]], matches_dict[dihedral[1]], matches_dict[dihedral[2]], matches_dict[dihedral[3]]
        dihedral[0], dihedral[1], dihedral[2], dihedral[3] = idx0, idx1, idx2, idx3
        dihedral[8] = basename(args.output) + "  " + str(old_idx0) + ' ' + str(old_idx1) + ' ' + str(old_idx2) + ' ' + str(old_idx3)
    itp.header.moleculetype.dihedrals.set_data(sorted(itp.header.moleculetype.dihedrals.data, key=lambda data: data[1]))

    # [exclusions]
    for exclusion in itp.header.moleculetype.exclusions.data:
        idx0, idx1 = matches_dict[exclusion[0]], matches_dict[exclusion[1]]
        exclusion[0] = min(idx0, idx1)
        exclusion[1] = max(idx0, idx1)
    itp.header.moleculetype.exclusions.data.sort()

    # Write out
    itp.write(args.output)

if __name__ == '__main__':
    main()
