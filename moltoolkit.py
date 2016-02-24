#!/usr/bin/env python3
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
    def __init__(self, mol_file):
        if not os.path.exists(mol_file):
            print('PDB file %s is not exist!' % mol_file)
            exit(1)
        self.__fileformat = mol_file[-3:]
        self.__filename = mol_file

        self.__obConversion = ob.OBConversion()
        self.__obConversion.SetInAndOutFormats(self.__fileformat, 'smi')
        self.__obmol = ob.OBMol()
        self.__obConversion.ReadFile(self.__obmol, self.__filename)

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
        self.most_coordinates
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
    def symmetrys(self):
        '''Find the atom symmetrys in the molecule.

        There are two versions to calculate the symmetrys of the molecule: C++ or Python.
        The former is more accurate. The latter will crash when there are two many atoms in the molecule.
        So the program will check if there is a C++ version, if not, using python.

        @Returns
            a tuple such as [[1,2,3], [4,5]]
        '''
        if not shutil.which('Symmetry'):
            return self.__get_symmetrys_python()
        else:
            return self.__get_symmetrys_cpp()

    def __get_symmetrys_python(self):
        '''Find the atom symmetrys in the molecule (Python Version).

        There are two versions to calculate the symmetrys of the molecule: C++ or Python.
        The former is more accurate. The latter will crash when there are two many atoms in the molecule.
        So the program will check if there is a C++ version, if not, using python.

        @Returns
            a tuple such as [(1,2,3), (4,5)]
        '''
        # Get symmetry classes
        gs = ob.OBGraphSym(self.__obmol)
        symclasses = ob.vectorUnsignedInt()
        gs.GetSymmetry(symclasses)

        # Get automorphisms
        automorphs = ob.vvpairUIntUInt()
        print("\nGet symmetrys python version: When there are too many atoms. The program will crash  as a error: memory exceed or timeout. Please don't use this result. Please use the cpp version.\n")
        ob.FindAutomorphisms(self.__obmol, automorphs, symclasses, ob.OBBitVec(), 40960000)
        amlist = []

        for x in automorphs:
            # openbabel有時原子索引從0開始，需要校驗一下
            for (y, z) in x:
                if y != z:
                    pair = sorted([y + 1, z + 1])
                    if not pair in amlist:
                        amlist.append(pair)
        amlist.sort()

        # Merge all the symmetrys into one list合并所有的對稱原子至一個小列表內
        last_amlist = []
        while len(amlist) > 0:
            same_sym = self.__find_same_item_list(amlist[0], amlist)
            temp_list = []
            for i in same_sym:
                temp_list += list(i)
                amlist.remove(list(i))
            last_amlist.append(sorted(list(set(temp_list))))
        return last_amlist

    def __get_symmetrys_cpp(self):
        '''Find the atom symmetrys in the molecule (C++ Version).

        There are two versions to calculate the symmetrys of the molecule: C++ or Python.
        The former is more accurate. The latter will crash when there are two many atoms in the molecule.
        So the program will check if there is a C++ version, if not, using python.

        It needs a executable program in path (Symmetry).
        @Returns
            a tuple such as [(1,2,3), (4,5)]
        '''
        result = subprocess.check_output(["Symmetry", self.__fileformat, self.__filename], universal_newlines=True)
        result1 = '[' + result.split('\n')[0] + ']'
        result_list = eval(result1)
        return result_list

    def __find_same_item_list(self, list_item, list_list):
        ''' 根据一个列表a，此列表a内包含有一个列表b，根据所给的参数列表c，循环c中的元素d，在b中查找是否与d相同的元素，如果有，b中的其它元素也要继续像a一样在所有小列表b中查找，就像一个多叉树一样#

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


    def getSubstructureMaps(self, mol2, includeAtoms=[]):
        '''Find the indice maps of substructure and the molecule.

        Args:
            mol2: The substructure of the mol.

        Returns:
            a matches list contains indices of the substructure and the moleculesuch such as [(1,4), (2, 5)], if none, return a empty list
        '''
        query = ob.CompileMoleculeQuery(mol2.__obmol)
        mapper = ob.OBIsomorphismMapper.GetInstance(query);
        if len(includeAtoms) == 0:
            maps = ob.vpairUIntUInt()
            mapper.MapFirst(self.__obmol, maps);
            return list([(i + 1, j + 1) for (i, j) in maps])
        else:
            maps = ob.vvpairUIntUInt()
            mapper.MapAll(self.__obmol, maps, ob.OBBitVec(), 40960000)
            for group in maps:
                # flag for if include_atoms is in the group
                notThis = False
                mol1List = [j for (i, j) in group]
                for i in includeAtoms:
                    if not i in mol1List:
                        notThis = True
                        break
                if notThis:
                   continue
                else:
                    return list([(i + 1, j + 1) for (i, j) in group])
            else:
                return []
        #for pair in matches:
        #    if mol2.get_atom(pair[0]).atomicNum != self.get_atom(pair[1]).atomicNum:
        #        return []
        #else:
        #    return matches

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
        print('分子對稱原子對:', mol.symmetrys)

    if args.identical:
        print('Identical Molecule Matches: ', mol.getIdenticalMaps(Mol(args.identical)))

    if args.substructuremaps:
        print('Substructure Matches: ', mol.getSubstructureMaps(Mol(args.substructuremaps)))

if __name__ == '__main__':
    main()
