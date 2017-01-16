#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Copyright (C) Hugh Gao
#
# This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along with this program. If not, see http://www.gnu.org/licenses/.
##
# @file moltoolkit.py
# @brief A toolkit about molecule.
# @author Hugh Gao, hughgao01@gmail.com
# @version 0.1
# @date 2012-04-21


import openbabel as ob
import os
import shutil
import subprocess
import argparse


class Mol:
    '''A toolkit to get the information of a molecule, manipulate atom and so on.'''
    def __init__(self, file_or_str, format="pdb"):
        self.__fileformat = format
        self.__filename = file_or_str
        self.__obConversion = ob.OBConversion()
        self.__obConversion.SetInAndOutFormats(self.__fileformat, 'smi')
        self.__obmol = ob.OBMol()
        if format == "pdb":
            if not os.path.exists(file_or_str):
                print('PDB file %s is not exist!' % file_or_str)
                exit(1)
            else:
                self.__obConversion.ReadFile(self.__obmol, self.__filename)
        elif format == "smi":
            self.__obConversion.ReadString(self.__obmol, self.__filename)

    @property
    def smiles(self):
        self.__obConversion.SetInAndOutFormats(self.__fileformat, 'smi')
        return self.__obConversion.WriteString(self.__obmol)

    @property
    def most_coordinates(self):
        '''3D most coordinates of the molecule.

        Returns:
            The list as following: (smallest X, smallest Y, smallest Z, largest X, largest Y, largest Z)
        '''
        try:
            if self.__smallestX:
                return (self.__smallestX, self.__smallestY, self.__smallestZ, self.__largestX, self.__largestY, self.__largestZ)
        except AttributeError:
            self.__smallestX = 10000000
            self.__smallestY = 10000000
            self.__smallestZ = 10000000
            self.__largestX = -10000000
            self.__largestY = -10000000
            self.__largestZ = -10000000
            for atom in self.atoms:
                x, y, z = atom.coordinates
                if x < self.__smallestX:
                    self.__smallestX = x
                if y < self.__smallestY:
                    self.__smallestY = y
                if z < self.__smallestZ:
                    self.__smallestZ = z
                if x > self.__largestX:
                    self.__largestX = x
                if y > self.__largestY:
                    self.__largestY = y
                if z > self.__largestZ:
                    self.__largestZ = z
        return (self.__smallestX, self.__smallestY, self.__smallestZ, self.__largestX, self.__largestY, self.__largestZ)

    @property
    def lenths_at_axis(self):
        '''Get the 3 coordinates which define space of the molecule occupled.

        Returns:
            a tuple as (x_dis, y_dis, z_dis)
        '''
        return (self.__largestX - self.__smallestX, self.__largestY - self.__smallestY,
                self.__largestZ - self.__smallestZ)

    def __str__(self):
        return 'Atoms Num: ' + str(self.num_atoms) + '\nmost coordinates: ' + str(self.most_coordinates) +  '\nVolumn: ' + str(self.lenths_at_axis)

    @property
    def num_atoms(self):
        '''The number of the atoms in the molecule.

        Returns:
            the number of the atoms of the molecule.
        '''
        return self.__obmol.NumAtoms()

    @property
    def atoms(self):
        '''The the atoms in the molecule.

        Returns:
            The list contains all the atoms in the molecule.
        '''
        return [Atom(i) for i in ob.OBMolAtomIter(self.__obmol)]

    def get_atom(self, idx):
        '''Get the atom by idx of the atom.

        Args:
            idx: The index of the atom.
            Warning: the idx starts from 1.

        Returns:
            atom
        '''
        try:
            return self.atoms[idx - 1]
        except IndexError:
            print("Error: the idx you give for the atom is out of range. idx=", idx)
            exit(1)

    def delete_atom(self, atom):
        '''Delete the atom.

        Args:
            atom: The atom to be removed.
        '''
        return self.__obmol.DeleteAtom(atom.atom)

    def get_atom_idx(self, atom):
        '''Get the idx of the atom.

        Args:
            atom: The atom.

        Returns:
            index of the atom
            Warning: the idx of the atom starts from 0.
        '''
        return self.atoms.index(atom)

    def get_distance(self, atom1, atom2):
        '''Get the distance between two atoms.

        Args:
            atom1: The first atom.
            atom2: The second atom.

        Returns:
            The distance between the two atoms.
        '''
        dises = map(lambda x: pow(x, 2), self.get_distance_at_axis(atom1, atom2))
        return pow(sum(dises), 0.5)

    def get_distance_at_axis(self, atom1, atom2):
        '''Get the distance between two atoms at axis direction.

        Args:
            atom1: The first atom.
            atom2: The second atom.

        Returns:
            A tuple as the distance between the two atoms at axis: (x_dis, y_dis, z_dis)
        '''
        return [abs(atom1.coordinates[i] - atom2.coordinates[i]) for i in range(3)]

    @property
    def symmetries(self):
        '''Find the atom symmetrys in the molecule.

        @Returns
            a tuple such as [[1,2,3], [4,5]]
        '''
        # Get symmetry classes
        gs = ob.OBGraphSym(self.__obmol)
        sym_classes = ob.vectorUnsignedInt()
        gs.GetSymmetry(sym_classes)
        sym_list = [x for x in sym_classes]

        sym_maps = {}
        i = 0
        for atom in ob.OBMolAtomIter(self.__obmol):
            sym_grp = sym_list[i]
            if sym_grp in sym_maps:
                sym_maps[sym_grp].append(atom.GetIdx())
            else:
                sym_maps[sym_grp] = [atom.GetIdx()]
            i += 1

        return [v for k, v in sym_maps.items()]


    def __find_same_item_list(self, list_item, list_list):
        ''' 根据一个列表a，此列表a内包含有一个列表b，根据所给的参数列表c，循环c中的元素d，在b中查找是否与d相同的元素，如果有，
        b中的其它元素也要继续像a一样在所有小列表b中查找，就像一个多叉树一样#

        Args
            list_item 要查找的元素所组成的列表
            list_list 嵌套列表

        Returns
            具有這個元素的所有小tuple
        '''
        same_list = []
        items = []
        small_li = list_item[:]
        big_li = list_list[:]
        # 循环查找small_li内的元素，直到其为空
        while len(small_li) > 0:
            # 把已经找过的元素放到一个列表内，一会儿在新small_li内把它们清除
            items.append(small_li[0])
            # 开始递归大列表
            for j in big_li:
                #items[-1]就是要找的元素
                if items[-1] in j:
                    # 储存找到的小列表
                    same_list.append(tuple(j))
                    # 把新元素放到要查找的列表内
                    small_li += j
                    small_li = list(set(small_li))
            # 去除已经找过的元素
            list_set = set(small_li)
            list_set.difference_update(set(items))
            small_li = list(list_set)
        return list(set(same_list))

    def getIdenticalMaps(self, mol2):
        '''Get the indices mathes of two identical molecules.

        Args:
            mol: The molecule of the second molecule.

        Returns:
            [()], 索引對應list
        '''
        return self.getSubstructureMaps(mol2)


    def getSubstructureMaps(self, mol2):
        '''Find the indice maps of substructure and the molecule.

        Args:
            mol2: The substructure of the mol.

        Returns:
            a matches list contains indices of the substructure and the moleculesuch such as [(1,4), (2, 5)], if none, return a empty list
        '''
        query = ob.CompileMoleculeQuery(mol2.__obmol)
        mapper = ob.OBIsomorphismMapper.GetInstance(query)
        # it can't be implement now
        # masks = ob.OBBitVec()
        # for a in includeAtoms:
        #    masks.SetBitOn(a)
        maps = ob.vpairUIntUInt()
        # mapper.MapFirst(self.__obmol, maps, masks)
        mapper.MapFirst(self.__obmol, maps)
        result = list([(i + 1, j + 1) for (i, j) in maps])
        # remove the ones which not is the same type atom
        result1 = [pair for pair in result if self.get_atom(pair[1]).type == mol2.get_atom(pair[0]).type]
        if len(result) != len(result1):
            print("Warining: because the substructure may have extra hydrogen atoms at the end, "
                  "so the matching atoms with different types was removed.")
        return result1

class Atom:
    def __init__(self, atom):
        '''
        Args:
            atom: the atom for OBAtom
        '''
        self.__atom = atom

    @property
    def obatom(self):
        '''The obatom instance for the atom.'''
        return self.__atom

    @property
    def type(self):
        '''Get atom type by index of the atom.

        Returns:
            The type string of the atom.
        '''
        return self.__atom.GetType()

    @property
    def idx(self):
        '''Get the idx of the atom in the molecule.

        Returns:
            The idx of the atom.
            Warning: the idx starts from 1.
        '''
        return self.__atom.GetIdx()

    @property
    def atomID(self):
        '''Get the atomID of the atom.

        The atom ID is often the third column of the pdb.
        '''
        return self.__atom.GetResidue().GetAtomID(self.__atom)

    @property
    def neighbours(self):
        '''Get the neighbour atoms of the atom.

        Returns:
            The list contains the neighbour atoms of the atom.
        '''
        return [Atom(i) for i in ob.OBAtomAtomIter(self.__atom)]

    @property
    def coordinates(self):
        '''Get the 3d coordinate of the atom.

        Returns:
            A tuple as (x, y, z).
        '''
        x = self.__atom.x()
        y = self.__atom.y()
        z = self.__atom.z()
        return (x, y, z)

    @property
    def atomicNum(self):
        '''Get the 3d coordinate of the atom.

        Returns:
            A tuple as (x, y, z).
        '''
        return self.__atom.GetAtomicNum()

    def is_hydrogen(self):
        '''If the atom is hydrogen atom.'''
        return self.__atom.IsHydrogen()


def main():
    parser = argparse.ArgumentParser(description='Get the volumn of a molecule, the coordinate of atom, the distance between two atoms, the axial distance at x, y, z between two atoms, the smallest and largest coordinate at x, y, z axel of a molecule')
    parser.add_argument('-v', '--volumn', action='store_true', help='Get the squre volumn of a molecule')
    parser.add_argument('-c', '--coordinate', type=int, help='Get the coordinate of a atom. You must provide the index of a atom')
    parser.add_argument('-d', '--distance', nargs=2, type=int, help='Get the distance between atoms, it must provide indices of two atoms')
    parser.add_argument('-a', '--axle_distance', nargs=2, type=int, help='Get the distance at three axles between atoms, it must provide indices of two atoms')
    parser.add_argument('-p', '--pole', action='store_true', help='Get the 6 values which express the smallest and largest coordinates of a molecule')
    parser.add_argument('--symmetry', action='store_true', help='Get the symmetrys of molecule')
    parser.add_argument('--identical', help='Get the matches of two identical molecules. Please give a molecule file as argument.')
    parser.add_argument('--substructuremaps', help='Get the indices matches of a molecule and its substructure. Please give a molecule file as argument.')
    parser.add_argument('molfile', help='molecule file')
    args = parser.parse_args()
    mol = Mol(args.molfile)
    # 默认输出
    print(mol)
    if args.coordinate:
        print('Coordinate of atom which index is %s: %s' % (args.coordinate, mol.get_atom(args.coordinate).coordinates))
    if args.distance:
        atom1 = mol.get_atom(args.distance[0])
        atom2 = mol.get_atom(args.distance[1])
        print('Distance between atoms %s and %s:' % (args.distance[0], args.distance[1]), mol.get_distance(atom1, atom2))
    if args.axle_distance:
        atom1 = mol.get_atom(args.axle_distance[0])
        atom2 = mol.get_atom(args.axle_distance[1])
        print('Distance at axles between atoms %s and %s:' % (args.axle_distance[0], args.axle_distance[1]), mol.get_distance_at_axis(atom1, atom2))
    if args.pole:
        print('Six extreme coordinate of molecule:', mol.most_coordinates)

    if args.symmetry:
        print('分子對稱原子對:', mol.symmetries)

    if args.identical:
        print('Identical Molecule Matches: ', mol.getIdenticalMaps(Mol(args.identical)))

    if args.substructuremaps:
        print('Substructure Matches: ', mol.getSubstructureMaps(Mol(args.substructuremaps)))

if __name__ == '__main__':
    main()
