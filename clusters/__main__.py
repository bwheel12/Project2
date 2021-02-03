import algs

with open('ligand_information.csv') as f:
    ligands_txt = f.readlines()
    
all_ligands = []
for i in range(1,len(ligands_txt)):
    print(i)
    all_ligands.append(algs.ligand(ligands_txt[i]))