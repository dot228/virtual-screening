from rdkit import Chem
from rdkit import DataStructs

from chemistry.reactions import reaction
from chemistry.filters import pass_basic_filters
from chemistry.descriptors import calc_admet
from chemistry.fingerprints import (
    morgan_fp,
    topo_fp
)

from config import (
    LEAD_SMILES,
    MIN_MORGAN_SIM,
    MIN_TOPO_SIM
)

lead_mol = Chem.MolFromSmiles(
    LEAD_SMILES
)

lead_morgan = morgan_fp(
    lead_mol
)

lead_topo = topo_fp(
    lead_mol
)

def process_pair(pair):

    acid_smiles, amine_smiles = pair

    acid = Chem.MolFromSmiles(acid_smiles)
    amine = Chem.MolFromSmiles(amine_smiles)

    if acid is None or amine is None:
        return []
    acid_smiles = Chem.MolToSmiles(acid, canonical=True)
    amine_smiles = Chem.MolToSmiles(amine, canonical=True)

    acid = Chem.MolFromSmiles(acid_smiles)
    amine = Chem.MolFromSmiles(amine_smiles)

    rows = []

    try:

        products = reaction.RunReactants(
            (acid, amine)
        )

        for product in products:

            mol = product[0]

            try:
                Chem.SanitizeMol(mol)
            except:
                continue

            if not pass_basic_filters(mol):
                continue

            fp_morgan = morgan_fp(mol)
            fp_topo = topo_fp(mol)

            sim_morgan = DataStructs.TanimotoSimilarity(
                lead_morgan,
                fp_morgan
            )

            sim_topo = DataStructs.FingerprintSimilarity(
                lead_topo,
                fp_topo
            )

            if sim_morgan < MIN_MORGAN_SIM:
                continue

            if sim_topo < MIN_TOPO_SIM:
                continue

            rows.append({
                "Lead": LEAD_SMILES,
                "Acid": acid_smiles,
                "Amine": amine_smiles,
                "GeneratedSMILES": Chem.MolToSmiles(mol, canonical=True),
                "MorganSimilarity": round(sim_morgan, 4),
                "TopologicalSimilarity": round(sim_topo, 4),
                **calc_admet(mol)
            })

    except Exception:
        pass

    return rows