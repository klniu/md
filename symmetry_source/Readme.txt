Symmetry is a program to caculate the symmetrys of a molecule. It depends on the openbabel library.

Example:

For methane whose SMILES is C, and its atoms index is [1,2,3,4]. 1 is C; 2,3,4 is H. So 2,3,4 are symmetric. Using this program will get [[1],[2,3,4]].

== usuage ==
Command:
<pre>
Symmetry MoleculeFileFormat InputFileName
</pre>

MoleculeFileFormat is the format of the molecule file, e.g. pdb, mol2

InputFileName is the file path of the molecule.

e.g. Symmetry pdb /home/user/methane.pdb
