import unittest

from moltoolkit import Mol

class TestMoltoolkit(unittest.TestCase):

    def test_symmetries(self):
        mol = Mol("CCCC","smi")
        self.assertListEqual(mol.symmetries, [[1,4], [2,3]])

        mol = Mol("test_resources/4c-16.pdb")
        symmetries = mol.symmetries
        self.assertEqual(len(symmetries), 40)
        self.assertListEqual(symmetries, [[1], [2], [3], [4], [5], [6], [7], [8], [9], [10], [11], [12], [13], [14],
                                          [15, 23], [16, 22], [17], [18], [19, 20, 21], [24], [25], [26], [27, 28, 29],
                                          [30, 31], [32, 33], [34, 35], [36, 37], [38, 39], [40, 41], [42, 43],
                                          [44, 45], [46, 47], [48, 49], [50, 51], [52], [53, 56], [54, 55], [57, 58],
                                          [59, 60], [61, 62, 63]])

    def test_getSubstructureMaps(self):
        mol = Mol("CCCC","smi")
        fragment = Mol("C", "smi")

        matches = mol.getSubstructureMaps(fragment)
        self.assertEqual(len(matches), 1)
        self.assertListEqual(matches, [(1, 4)])

        # mask
        # matches = mol.getSubstructureMaps(fragment, [2])
        # self.assertEqual(len(matches), 1)
        # self.assertListEqual(matches, [(1, 2)])

        # complex examples
        mol = Mol("test_resources/shuang.pdb")
        fragment = Mol("test_resources/shuang-fragment.pdb")
        matches = mol.getSubstructureMaps(fragment)
        self.assertEqual(len(matches), 36)
        self.assertListEqual(matches, [(1, 143), (2, 62), (4, 142), (5, 61), (6, 140), (7, 141), (8, 60),
                                       (9, 138), (10, 139), (11, 59), (12, 58), (13, 136), (14, 137), (15, 57),
                                       (16, 134), (17, 135), (18, 56), (19, 55), (20, 132), (21, 133), (22, 54),
                                       (23, 130), (24, 131), (25, 53), (26, 52), (27, 128), (28, 129), (29, 51),
                                       (30, 126), (31, 127), (32, 50), (33, 48), (34, 49), (35, 47), (37, 124),
                                       (38, 125)])


if __name__ == '__main__':
    unittest.main()

