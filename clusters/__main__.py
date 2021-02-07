from clusters import algs
import numpy
import pickle
from datetime import datetime

with open('ligand_information.csv') as f:
    ligands_txt = f.readlines()
    
all_ligands = []
for i in range(1,len(ligands_txt)):
    if i % 3 == 0:
        all_ligands.append(algs.ligand(ligands_txt[i]))

    
legit_cluster = algs.HierarchicalClustering('single-linkage',1)
legit_cluster.get_data(all_ligands)
legit_cluster.cluster()

with open("legit_hier_clust.obj",'wb') as cluster_file:
    pickle.dump(legit_cluster, cluster_file)
    cluster_file.close()

del legit_cluster

with open('legit_hier_clust.obj','rb') as f:
    unpickler = pickle.Unpickler(f)
    legit_cluster = unpickler.load()
    f.close()