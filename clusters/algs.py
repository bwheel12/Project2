import numpy
from datetime import datetime
import numba
import random


class HierarchicalClustering():
    def __init__(self,distance_paradigm,cluster_number):
        self.dist_prdgm = distance_paradigm
        self.dendogram  = [] #somethings a matric maybe?
        self.cluster_number = cluster_number

        
    def update_cluster(self,cluster_number):
        self.cluster_number = cluster_number
        
    def get_data(self,objects_to_cluster):
        self.to_cluster = objects_to_cluster
        self.data_expanded = numpy.zeros((len(self.to_cluster),1024))
        i = 0
        for x in self.to_cluster:
            self.data_expanded[i,x.OnBits] = 1
            i = i + 1
    
    @numba.jit(parallel = True)    
    def calculate_distances(self):
        num_data = int(len(self.to_cluster))
        print(num_data)
        self.all_distances = numpy.full((num_data,num_data),100)
        for i in numba.prange(num_data):
            for j in numba.prange(i+1,num_data):
                self.all_distances[i,j] = self.distance(i,j)
                #print(i, j)
                if i%100 == 0:
                    print("i is ",i,end = '\r')
        
      
    
    ##function to determine distance between clusters for single linkage returns the indices of the current dendogram level to be lumped together
    ## *****NOTE ***** this method is super slow. I should instead calculate all this distances once and then look them up according to the rules in the following three functions
    @numba.jit(parallel=True)
    def single_linkage(self,dendogram_level):
        unique_cluster_list = numpy.unique(self.dendogram[dendogram_level,:]) #get the list of unique clusters at this level of the dendogram
        #print(unique_cluster_list)
        num_unique_clusters = len(unique_cluster_list) #get the number of unique clusters at this level of the dendogram
        #print(num_unique_clusters)
        

        
        for i in numba.prange(num_unique_clusters):
            for j in numba.prange(i+1,num_unique_clusters):
                entries_i = numpy.where(self.dendogram[dendogram_level,:] == unique_cluster_list[i])[0]
                entries_j = numpy.where(self.dendogram[dendogram_level,:] == unique_cluster_list[j])[0]
                #num_i = len(entries_i)
                #num_j = len(entries_j)

                #print(num_i)
                #print(num_j)
                single_distance = 100
                
                for z in entries_i:
                    for zed in entries_j:
                        single_temp_distance = self.all_distances[z,zed]
                        if single_temp_distance < single_distance:
                             single_distance = single_temp_distance
                                
                self.temp_distances[i,j] = single_distance
                
                
               
                #print(short_temp_distances)
                #self.temp_distances[i,j,dendogram_level] = numpy.amin(short_temp_distances)

               
        #print(temp_distances)      
        #print(numpy.amin(temp_distances))
        smallest_temp_distance = numpy.amin(self.temp_distances)
        first_cluster_index = numpy.where(self.temp_distances == smallest_temp_distance)[0][0]
        second_cluster_index = numpy.where(self.temp_distances == smallest_temp_distance)[1][0]
        
        
        return numpy.hstack((numpy.where(self.dendogram[dendogram_level,:] == unique_cluster_list[first_cluster_index])[0],numpy.where(self.dendogram[dendogram_level,:] == unique_cluster_list[second_cluster_index])[0]))       
                
        
        
    #function to determine distance between clusters for complete linkage
    def complete_linkage():
        pass
        
    #function to determine distance between clusters for average linkage
    def average_linkage():
        pass
    
    ##function to look up the distance between two ligands from existing table OBSOLETE
    def distance_lookup(self,object_1_index,object_2_index):
        return self.all_distances[self.to_cluster.index(object_1),self.to_cluster.index(object_2)]
    
    
    #function to determine distance between two ligands, uses euclidian distance, accepts ligand objects
    def distance(self,object_1_index,object_2_index):
        total_distance = len(self.to_cluster[object_1_index].OnBits) + len(self.to_cluster[object_2_index].OnBits) -  2*numpy.sum(self.to_cluster[object_1_index].OnBits == self.to_cluster[object_2_index].OnBits)
        euclid_dist = numpy.sqrt(total_distance)
        #print(euclid_dist)
        return euclid_dist
        
    ##function to determine clusters
    def cluster(self):
        self.calculate_distances()
        print('Distances_Done')
        self.dendogram = numpy.zeros((len(self.to_cluster),len(self.to_cluster)))
        self.dendogram[0,:] = numpy.arange(len(self.to_cluster)) #every object gets its own cluster number
        
        
        ## goes through 
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
    ##this class will implement traditional K-means clustering
    def __init__(self,cluster_number,OnBit_Dim):
        self.cluster_number = cluster_number
        self.Bit_Num        = OnBit_Dim

    def get_data(self,objects_to_cluster):
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
        self.similarity_matrix = numpy.zeros((self.total_objects,self.total_objects))
        for i in range(self.total_objects):
            for j in range(self.total_objects):
                self.similarity_matrix[i,j] = self.tanamota_coeff(self.to_cluster[i],self.to_cluster[j])
            print('i is',i,end='\r')
                
    def tanamota_coeff(self,object_1,object_2):
        intersection = numpy.sum(object_1.OnBits == object_2.OnBits)
        union = len(object_1.OnBits) + len(object_2.OnBits) - numpy.sum(object_1.OnBits == object_2.OnBits)
        return intersection/union
        
        
    #function to determine distance between two ligands, uses euclidian distance, expects 1D arrays of numerical values
    def distance(self,object_1,object_2):
        total_distance = 0
        obj_dimension = object_1.shape[0]
        for i in range(obj_dimension):
            total_distance = total_distance + (object_1[i]-object_2[i])**2
        euclid_dist = numpy.sqrt(total_distance)
        return euclid_dist
    
    def initialize_centroids(self):
        self.centroids = numpy.zeros((self.cluster_number,self.total_objects))
        available_indices = list(range(self.total_objects))
        for i in range(self.cluster_number):
            rand_location = random.randrange(len(available_indices))
            self.centroids[i,:] = self.similarity_matrix[available_indices[rand_location],:]
            available_indices.remove(available_indices[rand_location])
            #print(available_indices)
        #print(self.centroids)    
    
    def update_centroids(self):
        for i in range(self.cluster_number):
            cluster_data_expanded = self.similarity_matrix[numpy.where(self.cluster_assignments == i)[0],:]
            self.centroids[i,:] = cluster_data_expanded.mean(axis=0)
            if numpy.isnan(self.centroids[i,:]).any():
                rand_location = random.randrange(self.total_objects)
                self.centroids[i,:] = self.similarity_matrix[rand_location,:]
                print("Na fixed")
    
    def assign_cluster(self):
        temp_center_distances = numpy.zeros((self.total_objects,self.cluster_number))
        for i in range(self.total_objects):
            for j in range(self.cluster_number):
                temp_center_distances[i,j] = self.distance(self.similarity_matrix[i,:],self.centroids[j,:])
            self.cluster_assignments[i] = numpy.argmin(temp_center_distances[i,:])
        unique, counts = numpy.unique(self.cluster_assignments, return_counts=True)
        print(dict(zip(unique, counts)))
        #print(temp_center_distances)

    
    def is_cluster_over(self):
        
        if numpy.array_equal(self.centroids,self.centroids_1back) or numpy.array_equal(self.centroids,self.centroids_2back):
            self.keep_clustering = False
        
    
    def cluster(self):
        self.initialize_centroids()
        self.keep_clustering = True
        self.centroids_1back = numpy.zeros((self.cluster_number,self.total_objects))
        self.centroids_2back = numpy.zeros((self.cluster_number,self.total_objects))
        i = 1
        while self.keep_clustering:
            self.assign_cluster()
            self.centroids_2back = numpy.copy(self.centroids_1back)
            #print(self.centroids_2back)
            self.centroids_1back = numpy.copy(self.centroids)
            #print(self.centroids_1back)
            self.update_centroids()
            #print(self.centroids)
            
            print(self.distance(self.centroids[1,:],self.centroids_1back[1,:]))
            print(self.distance(self.centroids[2,:],self.centroids_1back[2,:]))
            print(self.distance(self.centroids[3,:],self.centroids_1back[3,:]))
            
            #self.keep_clustering = False
            self.is_cluster_over() 
            print('Iterrations Completed:',i,end='\r')
            i = i +1
        

    

class ligand():
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
      
        #print(self.ID,self.score,self.smiles,self.OnBits,end='\r')
