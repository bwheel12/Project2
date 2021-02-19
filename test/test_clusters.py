from clusters import algs
import numpy
import pytest

@pytest.fixture
def ligand_test():
    """This will read in the first 10 ligands from ligand_information.csv to use a dummy data set"""
    with open('ligand_information.csv') as f:
        ligands_txt = f.readlines()
    
    test_ligands = []
    for i in range(1,11):
        test_ligands.append(algs.ligand(ligands_txt[i]))
    
    return test_ligands

                   
                   
                   
def test_ligand_io(ligand_test):
    """this function test the creation of ligand objects"""
    assert ligand_test[0].ID == '0'
    assert ligand_test[0].score == -1.3
    assert ligand_test[0].smiles == 'N#C' 
    assert numpy.allclose(ligand_test[0].OnBits, numpy.array([360,489,915]))
                   
                   

def test_Hierarchical(ligand_test):
    """this function will test the k-means clustering for a series of k values and confirm that the end result has at least 1 ligand per cluster"""
    test_cluster = algs.HierarchicalClustering('single-linkage',1)
    test_cluster.get_data(ligand_test)
    test_cluster.cluster()
    for i in range(10):
        number_clusters = len(numpy.unique(test_cluster.dendogram[i,:]))
        assert number_clusters == (10-i)

def test_Partitioning(ligand_test):
    """this will cluster all 10 given objects and confirm that there are the appropriate number of clusters at each level of the dendogram"""               
    test_km_cluster = algs.PartitionClustering(10,1024)
    test_km_cluster.get_data(ligand_test)
    for i in range(1,10):
        test_km_cluster.update_cluster_number(i)
        test_km_cluster.cluster()
        number_clusters = len(numpy.unique(test_km_cluster.cluster_assignments))
        assert number_clusters == i


    
def test_Jaccard_Index(ligand_test):
    """This functions tests if identical clusters return a value of 1 and if completely different clusters return a value of 0"""
    test_km_cluster = algs.PartitionClustering(2,1024)
    test_km_cluster.get_data(ligand_test)
    test_km_cluster.cluster()
    wrong_clusters = numpy.array([20,21,22,23,24,25,26,27,28,39])
    correct_index = algs.Jaccard_Index(test_km_cluster.cluster_assignments,test_km_cluster.cluster_assignments)
    wrong_index   = algs.Jaccard_Index(test_km_cluster.cluster_assignments,wrong_clusters)
    assert correct_index == 1
    assert wrong_index   == 0
    
    
def test_silhouette_coeff(ligand_test):
    """This function tests if coeff falls in expected range"""
    test_km_cluster = algs.PartitionClustering(2,1024)
    test_km_cluster.get_data(ligand_test)
    test_km_cluster.cluster()
    distance_matrix_test = algs.euclid_distance(test_km_cluster.similarity_matrix)
    for j in range(10):
        test_score = algs.silhouette_score(distance_matrix_test,numpy.transpose(test_km_cluster.cluster_assignments)[0],j)
    
        assert test_score > -1 and test_score < 1