import unittest
import project4
import data_for_test

class Project4Test(unittest.TestCase):
    def test_protein_seq_creator(self):
        pdb_name = '7F5U'
        RE = project4.protein_seq_creator(pdb_name)
        TE = data_for_test.protein
        self.assertEqual(TE, RE)

    def test_surface_analysis(self):
        sequence = data_for_test.protein
        RE = project4.surface_analysis(sequence)
        TE = data_for_test.resav
        self.assertEqual(TE, RE)

    def test_site_seq_creator(self):
        prot_site = 'F/Y'
        RE = project4.site_seq_creator(prot_site)
        TE = 'FY', 'F', 'Y'
        self.assertEqual(TE, RE)

    def test_proteolysis_sites(self):
        seq = data_for_test.protein
        site_seq = 'FY'
        site_seq_left = 'F'
        site_seq_right = 'Y'
        RE = project4.proteolysis_sites(seq, site_seq, site_seq_left, site_seq_right)
        TE = {1: ['122/123', '807/808', '1492/1493', '2177/2178']}
        self.assertEqual(TE, RE)

if __name__ == '__main__':
    unittest.main()
