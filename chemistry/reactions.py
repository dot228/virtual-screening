from rdkit.Chem import AllChem

reaction = AllChem.ReactionFromSmarts(
    "[OX2H:1].[CH2:2][Cl,Br]>>[O:1][CH2:2]"
)