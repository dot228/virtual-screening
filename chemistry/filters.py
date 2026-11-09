from rdkit.Chem import Descriptors

def pass_basic_filters(mol):
    mw = Descriptors.MolWt(mol)
    if mw > 500:
        return False

    logp = Descriptors.MolLogP(mol)
    if logp > 5:
        return False

    return True