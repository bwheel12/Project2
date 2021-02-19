import numpy
from datetime import datetime
import random
from scipy.spatial import distance as scipy_dist

class HierarchicalClustering():
    """This class will create an object with which has the necessary attributes and functions to perform bottom-up agglomerative hierarchical clustering on a set of provided ligand objects."""
    def __init__(self,distance_paradigm,cluster_number):
        self.dist_prdgm = distance_paradigm
        self.dendogram  = [] 
        self.cluster_number = cluster_number

        
    def update_cluster(self,cluster_number):
        """This function resets the number of clusters to find in the algorithm by the number given by the user here"""
        self.cluster_number = cluster_number
        
    def get_data(self,objects_to_cluster):
        """This function takes the list of ligand objects to cluster given by the user and assigns them to the internal variable to_cluster.
            This function also uses the given ligands to populate the expanded data_expanded matrix"""
        self.to_cluster = objects_to_cluster
        self.data_expanded = numpy.zeros((len(self.to_cluster),1024))
        i = 0
        for x in self.to_cluster:
            self.data_expanded[i,x.OnBits] = 1
            i = i + 1
    
    def calculate_distances(self):
        """This function is called internally when the cluster function is called by the user. It calls the distance fuction for each ligand pair.
        It populates this data into the all_data internal attribute."""
        num_data = int(len(self.to_cluster))
        print(num_data)
        self.all_distances = numpy.full((num_data,num_data),100)
        for i in range(num_data):
            for j in range(i+1,num_data):
                self.all_distances[i,j] = self.distance(i,j)
                #print(i, j)
                if i%100 == 0:
                    print("i is ",i,end = '\r')
        
      
    
    def single_linkage(self,dendogram_level):
        """function to determine distance between clusters for single linkage returns the indices of the current dendogram level to be lumped together"""
        unique_cluster_list = numpy.unique(self.dendogram[dendogram_level,:]) #get the list of unique clusters at this level of the dendogram
        num_unique_clusters = len(unique_cluster_list) #get the number of unique clusters at this level of the dendogram
        

        
        for i in range(num_unique_clusters):
            for j in range(i+1,num_unique_clusters):
                entries_i = numpy.where(self.dendogram[dendogram_level,:] == unique_cluster_list[i])[0]
                entries_j = numpy.where(self.dendogram[dendogram_level,:] == unique_cluster_list[j])[0]

                single_distance = 100
                
                for z in entries_i:
                    for zed in entries_j:
                        single_temp_distance = self.all_distances[z,zed]
                        if single_temp_distance < single_distance:
                             single_distance = single_temp_distance
                                
                self.temp_distances[i,j] = single_distance
                
        smallest_temp_distance = numpy.amin(self.temp_distances)
        first_cluster_index = numpy.where(self.temp_distances == smallest_temp_distance)[0][0]
        second_cluster_index = numpy.where(self.temp_distances == smallest_temp_distance)[1][0]
        
        
        return numpy.hstack((numpy.where(self.dendogram[dendogram_level,:] == unique_cluster_list[first_cluster_index])[0],numpy.where(self.dendogram[dendogram_level,:] == unique_cluster_list[second_cluster_index])[0]))       
    
    
    def distance(self,object_1_index,object_2_index):
        """function to determine distance between two ligands, uses euclidian distance"""
        total_distance = len(self.to_cluster[object_1_index].OnBits) + len(self.to_cluster[object_2_index].OnBits) -  2*numpy.sum(self.to_cluster[object_1_index].OnBits == self.to_cluster[object_2_index].OnBits)
        euclid_dist = numpy.sqrt(total_distance)
        #print(euclid_dist)
        return euclid_dist
        
    def cluster(self):
        """function to determine clusters"""
        self.calculate_distances()
        print('Distances_Done')
        self.dendogram = numpy.zeros((len(self.to_cluster),len(self.to_cluster)))
        self.dendogram[0,:] = numpy.arange(len(self.to_cluster)) #every object gets its own cluster number
        
        
        ## goes through agglomeration procedure
        for zed in range(1,len(self.to_cluster)-self.cluster_number+1):
            self.temp_distances = numpy.full((len(self.to_cluster),len(self.to_cluster)),100) # square matrix to find distance between all existing clusters
            indices_to_lump = self.single_linkage(zed-1) #return an array of indices in the dendogram to lump together
            
            lowest_previous_cluster = numpy.amin(self.dendogram[zed-1,indices_to_lump])
            self.dendogram[zed,:] = self.dendogram[zed-1,:]
            self.dendogram[zed,indices_to_lump] = lowest_previous_cluster

            
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print("Current Time 2 =", current_time, "zed is ", zed,end='\r')

        
        
    
    
class PartitionClustering():
    """This class implements vanilla K-means clustering"""
    
    
    def __init__(self,cluster_number,OnBit_Dim):
        """Two parameters must be passed upon class initialization
        cluster_number: int, determines the number of clusters to find
        OnBit_Dim: int, indicate the dimension space of the OnBits so they can be appropriately unpacked
        """
        self.cluster_number = cluster_number
        self.Bit_Num        = OnBit_Dim

    def get_data(self,objects_to_cluster):
        """This function reads in the ligand objects and sets up essential class attributes to be used in clustering
        """
        self.to_cluster = objects_to_cluster
        self.total_objects = len(self.to_cluster)
        self.cluster_assignments = numpy.zeros((self.total_objects,1))
        self.data_expanded = numpy.zeros((self.total_objects,self.Bit_Num))
        i = 0
        for x in self.to_cluster:
            self.data_expanded[i,x.OnBits] = 1
            i = i + 1
        self.create_similarity_matrix()
        print('Data Imported')

            
    def create_similarity_matrix(self):
        """This function goes through each pair of ligands in to_cluster and calculates their similarity by the tanimoto coefficient and assigns this 
        to the internal similarity matrix attribute."""
        self.similarity_matrix = numpy.zeros((self.total_objects,self.total_objects))
        for i in range(self.total_objects):
            for j in range(self.total_objects):
                self.similarity_matrix[i,j] = self.tanamota_coeff(self.to_cluster[i],self.to_cluster[j])
            print('i is',i,end='\r')
                
    def tanamota_coeff(self,object_1,object_2):
        """This function determines the similarity of two ligand objects by calculating their tanimoto coefficient."""
        intersection = numpy.sum(object_1.OnBits == object_2.OnBits)
        union = len(object_1.OnBits) + len(object_2.OnBits) - numpy.sum(object_1.OnBits == object_2.OnBits)
        return intersection/union
        
        
    def distance(self,object_1,object_2):
        """function to determine distance between two ligands, uses euclidian distance, expects 1D arrays of numerical values"""
        total_distance = 0
        obj_dimension = object_1.shape[0]
        for i in range(obj_dimension):
            total_distance = total_distance + (object_1[i]-object_2[i])**2
        euclid_dist = numpy.sqrt(total_distance)
        return euclid_dist
    
    def initialize_centroids(self):
        """This function initializes the number of centroids indicated by cluster_number as random ligands in the similarity space"""
        self.centroids = numpy.zeros((self.cluster_number,self.total_objects))
        available_indices = list(range(self.total_objects)) #used to make sure centroids aren't initialized to the same ligand
        for i in range(self.cluster_number):
            rand_location = random.randrange(len(available_indices))
            self.centroids[i,:] = self.similarity_matrix[available_indices[rand_location],:]
            available_indices.remove(available_indices[rand_location])  
    
    def update_centroids(self):
        """This function finds the average location of ligands assigned to a given cluster and re-assigns the centroid to that location.
        If there is an emtpy cluster the centriod is reassigned to a random ligand in the data set."""
        for i in range(self.cluster_number):
            cluster_data_expanded = self.similarity_matrix[numpy.where(self.cluster_assignments == i)[0],:]
            self.centroids[i,:] = cluster_data_expanded.mean(axis=0)
            if numpy.isnan(self.centroids[i,:]).any():
                rand_location = random.randrange(self.total_objects)
                self.centroids[i,:] = self.similarity_matrix[rand_location,:]
                print("Na fixed")
    
    def assign_cluster(self):
        """This function finds the closest centroid to each ligand. It then assigns that centroid number as the cluster for that ligand"""
        temp_center_distances = numpy.zeros((self.total_objects,self.cluster_number))
        for i in range(self.total_objects):
            for j in range(self.cluster_number):
                temp_center_distances[i,j] = self.distance(self.similarity_matrix[i,:],self.centroids[j,:])
            self.cluster_assignments[i] = numpy.argmin(temp_center_distances[i,:])
        unique, counts = numpy.unique(self.cluster_assignments, return_counts=True)
        print(dict(zip(unique, counts)))
        

        
    def update_cluster_number(self,cluster_number):
        """This function updates the number of clusters to find"""
        self.cluster_number = cluster_number
    
    def is_cluster_over(self):
        """This function determines whether to stop the clustering procedure. If the current centroids are the same as either the previous
        or 2nd back set of centroids then clustering stops."""
        
        if numpy.array_equal(self.centroids,self.centroids_1back) or numpy.array_equal(self.centroids,self.centroids_2back):
            self.keep_clustering = False
        
    
    def cluster(self):
        """This function performs all the necessary steps to execute k-means clustering"""
        self.initialize_centroids()
        self.keep_clustering = True
        self.centroids_1back = numpy.zeros((self.cluster_number,self.total_objects))
        self.centroids_2back = numpy.zeros((self.cluster_number,self.total_objects))
        i = 1
        while self.keep_clustering:
            self.assign_cluster()
            self.centroids_2back = numpy.copy(self.centroids_1back)
            self.centroids_1back = numpy.copy(self.centroids)
            self.update_centroids()

            self.is_cluster_over() 
            print('Iterrations Completed:',i,end='\r')
            i = i +1
        

    

class ligand():
    """This simple class organizes all the import attributes for a ligand into a single object with each held as an attribute"""
    def __init__(self,ligand_string):
        import numpy
        first_comma = ligand_string.find(",")
        second_comma = ligand_string.find(",",first_comma+1)
        third_comma  = ligand_string.find(",",second_comma+1)
        final_n      = ligand_string.find("\n",second_comma+1)
        
        self.ID      = ligand_string[0:first_comma]
        self.score   = float(ligand_string[first_comma+1:second_comma])
        self.smiles  = ligand_string[second_comma+1:third_comma]
        onBits       = ligand_string[third_comma+2:final_n-1]
        self.OnBits  = numpy.array([int(i) for i in onBits.split(',')])
      

        
##below are functions that are useful for evaluating clusters

def silhouette_score(distance_matrix,cluster_list,POI):
    """this defines the silhouette score for a given point"""
    a_temp = 0
    b_temp = 0 
    score  = 0
    close_cluster = closest_cluster(distance_matrix,cluster_list,POI) #first find the closest cluster
    parent_cluster_indices = numpy.where(cluster_list == cluster_list[POI])[0] #get the ligand indices for the ligands in the same cluster
    close_cluster_indices   = numpy.where(cluster_list == close_cluster)[0] #get the ligand indices for the ligands in the closest cluster
    num_parent_cluster = len(parent_cluster_indices)
    num_close_cluster  = len(close_cluster_indices)
    
    #find the average distance to ligands in the same cluster
    for i in range(num_parent_cluster):
        if parent_cluster_indices[i] != POI:
            a_temp = a_temp + distance_matrix[parent_cluster_indices[i],POI]
    if num_parent_cluster > 1:
        a_final = a_temp/(num_parent_cluster-1)
    
    #if cluster is a single ligand define as 0
    else:
        a_final = 0
    
    #find the average distance to all the ligands in the closest cluster
    for j in range(num_close_cluster):
        b_temp = b_temp + distance_matrix[close_cluster_indices[j],POI]
    b_final = b_temp/num_close_cluster
    
    score = (b_final - a_final)/max(a_final,b_final)
    
    
    return score


def Jaccard_Index(cluster_list_1,cluster_list_2):
    """Compares the similarity in clustering results by calculating the jaccard index"""
    num_ligand = len(cluster_list_1)
    f_11 = 0 # same in both
    f_01 = 0 # same in 1 not in 2
    f_10 = 0 # same in 2 not in 1
    f_00 = 0 # diff in both
    
    for i in range(0,num_ligand):
        for j in range(i+1,num_ligand):
            ligand_1_clust_1 = cluster_list_1[i]
            ligand_1_clust_2 = cluster_list_2[i]
            ligand_2_clust_1 = cluster_list_1[j]
            ligand_2_clust_2 = cluster_list_2[j]
            
            if ligand_1_clust_1 == ligand_2_clust_1 and ligand_1_clust_2 == ligand_2_clust_2:
                f_11 = f_11 + 1
            if ligand_1_clust_1 == ligand_2_clust_1 and ligand_1_clust_2 != ligand_2_clust_2:
                f_01 = f_01 + 1
            if ligand_1_clust_1 != ligand_2_clust_1 and ligand_1_clust_2 == ligand_2_clust_2:
                f_10 = f_10 + 1
            if ligand_1_clust_1 != ligand_2_clust_1 and ligand_1_clust_2 != ligand_2_clust_2:
                f_00 = f_00 + 1
    
    print(f_11,f_01,f_10,f_00)
    jaccard_index = f_11/(f_01+f_10+f_11)
    return jaccard_index


def euclid_distance(similarity_matrix):
    """This function takes a similarity matrix and creates a euclidean distance matrix based upon that
    """
    obj_dimension = similarity_matrix.shape[0]
    print(obj_dimension)
    distance_matrix = numpy.zeros((obj_dimension,obj_dimension))
    for m in range(obj_dimension):
        for n in range(obj_dimension):
            distance_matrix[m,n] = scipy_dist.euclidean(similarity_matrix[m,:],similarity_matrix[n,:])
        print(m,end='\r')
    return distance_matrix


def closest_cluster(distance_matrix,cluster_list,ligand_index):
    """This function finds the cluster that is closest to the given ligand that is not contained within the ligands home cluster"""
    num_ligands = distance_matrix.shape[0]
    current_cluster = cluster_list[ligand_index]
    closest_distance = num_ligands #similarity scores are bounded at 1. So max euclidean distance is sqrt(num_ligands). Therefore, num_ligands linearly should always exceed 
    close_cluster = -1 #no cluster is ever called -1, so returning -1 indicates an error
    for i in range(num_ligands):
        temp_distance = distance_matrix[ligand_index,i]
        if temp_distance < closest_distance and cluster_list[i] != current_cluster:
            closest_distance = temp_distance
            close_cluster = cluster_list[i]
    return close_cluster