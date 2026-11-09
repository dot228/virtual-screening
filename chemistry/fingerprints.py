from rdkit import Chem
from rdkit.Chem import AllChem
def morgan_fp(mol):
    return AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=2048)
def topo_fp(mol):
    return Chem.RDKFingerprint(mol)