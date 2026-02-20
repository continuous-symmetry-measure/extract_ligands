# tests/test_extract_ligands.py
import pytest
from src.extract_ligands.extrac_ligands import is_het, is_main_ligand


class TestIsHet:
    """Test HETATM identification"""

    def test_is_het_for_hetatm(self):
        """HETATM should return True"""

        class MockResidue:
            id = ('H_XYZ', 100, ' ')

        assert is_het(MockResidue()) is True

    def test_is_het_for_water(self):
        """Water should return False"""

        class MockResidue:
            id = ('W', 100, ' ')

        assert is_het(MockResidue()) is False

    def test_is_het_for_standard_residue(self):
        """Standard amino acid should return False"""

        class MockResidue:
            id = (' ', 100, ' ')

        assert is_het(MockResidue()) is False


class TestIsMainLigand:
    """Test ligand filtering logic"""

    def test_main_ligand_not_in_filters(self):
        """Unknown ligand should be considered main ligand"""

        class MockResidue:
            def get_resname(self):
                return 'XYZ'

        assert is_main_ligand(MockResidue()) is True

    def test_ion_not_main_ligand(self):
        """Ions should not be main ligands"""

        class MockResidue:
            def get_resname(self):
                return 'CA'  # Calcium ion

        assert is_main_ligand(MockResidue()) is False

    def test_multiple_ions_not_main_ligand(self):
        """Test various ions are filtered"""
        ions_to_test = ['ZN', 'MG', 'FE', 'NA', 'K']

        for ion in ions_to_test:
            class MockResidue:
                def get_resname(self):
                    return ion

            assert is_main_ligand(MockResidue()) is False, f"{ion} should not be main ligand"

    def test_small_molecule_not_main_ligand(self):
        """Small molecules should not be main ligands"""

        class MockResidue:
            def get_resname(self):
                return 'EDO'  # Ethylene glycol

        assert is_main_ligand(MockResidue()) is False

    def test_multiple_small_mols_not_main_ligand(self):
        """Test various small molecules are filtered"""
        small_mols_to_test = ['GOL', 'SO4', 'PO4', 'ACT', 'DMS']

        for mol in small_mols_to_test:
            class MockResidue:
                def get_resname(self):
                    return mol

            assert is_main_ligand(MockResidue()) is False, f"{mol} should not be main ligand"

    def test_non_standard_residue_not_main_ligand(self):
        """Non-standard amino acids should not be main ligands"""

        class MockResidue:
            def get_resname(self):
                return 'MSE'  # Selenomethionine

        assert is_main_ligand(MockResidue()) is False

    def test_multiple_non_standard_res_not_main_ligand(self):
        """Test various non-standard residues are filtered"""
        non_std_to_test = ['MSE', 'CSO', 'AIB', 'NLE']

        for res in non_std_to_test:
            class MockResidue:
                def get_resname(self):
                    return res

            assert is_main_ligand(MockResidue()) is False, f"{res} should not be main ligand"

    def test_actual_drug_ligand_is_main(self):
        """Real drug-like ligands should be main ligands"""
        # These are typical 3-letter codes NOT in any filter list
        drug_ligands = ['ATP', 'HEM', 'FAD', 'NAD', 'GTP']

        for ligand in drug_ligands:
            class MockResidue:
                def get_resname(self):
                    return ligand

            assert is_main_ligand(MockResidue()) is True, f"{ligand} should be main ligand"

    def test_case_sensitive_ligand_name(self):
        """Ligand names are case-sensitive"""

        class MockResidue:
            def get_resname(self):
                return 'ca'  # lowercase - not in ions list (which has 'CA')

        # This tests if your filter is case-sensitive
        assert is_main_ligand(MockResidue()) is True


class TestIsMainLigandEdgeCases:
    """Test edge cases for ligand filtering"""

    def test_empty_residue_name(self):
        """Empty residue name should be main ligand"""

        class MockResidue:
            def get_resname(self):
                return ''

        assert is_main_ligand(MockResidue()) is True

    def test_longer_residue_name(self):
        """Longer residue names not in filters should be main ligands"""

        class MockResidue:
            def get_resname(self):
                return 'LONGNAME'

        assert is_main_ligand(MockResidue()) is True