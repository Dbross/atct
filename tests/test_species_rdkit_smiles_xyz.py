"""Species.to_rdkit_mol: SMILES/InChI topology + ATcT XYZ coordinates (no HTTP)."""

import unittest
from unittest import mock

try:
    from rdkit import Chem

    RDKIT = True
except ImportError:
    RDKIT = False

from atct.api.models import Species, XYZData


# Coordinates aligned with artifacts/ga_crossval/molfiles_atct_xyz/_117469-28-0_0.mol (B3LYP frame).
_CYCLPROPENYL_ANION_XYZ = """7
test
C         -1.587029      -0.012887       0.000000
C         -0.238429      -0.006279       0.000000
C          1.008971      -0.605301       0.000000
C          1.044927       0.753778       0.000000
H         -2.152722       0.915358       0.000000
H         -2.161938      -0.936393       0.000000
H          1.594668      -1.518634       0.000000
""".strip().splitlines()


@unittest.skipUnless(RDKIT, "RDKit not installed")
class TestSpeciesRdkitSmilesXyz(unittest.TestCase):
    def test_charged_cyclopropenyl_roundtrips_molblock(self) -> None:
        sp = Species(
            atct_tn_version=None,
            atct_id="*117469-28-0*0",
            name=None,
            formula="C4H3-",
            delta_h_0k=None,
            delta_h_298k=None,
            delta_h_298k_uncertainty=None,
            smiles="C=C1C=[C-]1",
            casrn="*117469-28-0",
            inchi="InChI=1S/C4H3/c1-4-2-3-4/h2H,1H2/q-1",
            inchi_key="OGLPWLKRCOJIET-UHFFFAOYSA-N",
            charge=-1,
            xyz=list(_CYCLPROPENYL_ANION_XYZ),
        )
        mol = sp.to_rdkit_mol()
        self.assertEqual(mol.GetNumAtoms(), 7)
        mb = Chem.MolToMolBlock(mol)
        back = Chem.MolFromMolBlock(mb, sanitize=True)
        self.assertIsNotNone(back)
        ik = Chem.MolToInchiKey(mol)
        self.assertEqual(ik, "OGLPWLKRCOJIET-UHFFFAOYSA-N")

    @mock.patch("rdkit.Chem.AllChem.EmbedMolecule", return_value=-1)
    def test_charged_cyclopropenyl_when_embed_fails_uses_2d_distance_template(self, _mock_embed) -> None:
        """If 3D embedding never succeeds, 2D coords still allow XYZ permutation."""
        sp = Species(
            atct_tn_version=None,
            atct_id="*117469-28-0*0",
            name=None,
            formula="C4H3-",
            delta_h_0k=None,
            delta_h_298k=None,
            delta_h_298k_uncertainty=None,
            smiles="C=C1C=[C-]1",
            casrn="*117469-28-0",
            inchi="InChI=1S/C4H3/c1-4-2-3-4/h2H,1H2/q-1",
            inchi_key="OGLPWLKRCOJIET-UHFFFAOYSA-N",
            charge=-1,
            xyz=list(_CYCLPROPENYL_ANION_XYZ),
        )
        mol = sp.to_rdkit_mol()
        self.assertEqual(mol.GetNumAtoms(), 7)
        mb = Chem.MolToMolBlock(mol)
        back = Chem.MolFromMolBlock(mb, sanitize=True)
        self.assertIsNotNone(back)
        self.assertEqual(Chem.MolToInchiKey(mol), "OGLPWLKRCOJIET-UHFFFAOYSA-N")


if __name__ == "__main__":
    unittest.main()
