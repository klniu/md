import unittest

from moltoolkit import Mol

class TestMoltoolkit(unittest.TestCase):

    def test_symmetries(self):
        mol = Mol("smi", "CCCC")
        self.assertListEqual(mol.symmetries, [[1,4], [2,3]])

        mol = Mol("pdb", "test_resources/4c-16.pdb")
        symmetries = mol.symmetries
        self.assertEqual(len(symmetries), 40)
        self.assertListEqual(symmetries, [[53, 56], [54, 55], [19, 20, 21], [44, 45], [42, 43], [40, 41], [38, 39],
                                          [30, 31], [34, 35], [36, 37], [27, 28, 29], [32, 33], [46, 47], [61, 62, 63],
                                          [50, 51], [57, 58], [48, 49], [59, 60], [52], [16, 22], [15, 23], [17], [14],
                                          [18], [8], [7], [1], [5], [3], [6], [4], [2], [9], [11], [25], [13], [10],
                                          [26], [12], [24]])


if __name__ == '__main__':
    unittest.main()

