import numpy
from datetime import datetime



class HierarchicalClustering():
    def __init__(self,distance_paradigm,cluster_number):
        self.dist_prdgm = distance_paradigm
        self.dendogram  = [] #somethings a matric maybe?
        self.cluster_number = cluster_number

        
    def update_cluster(self,cluster_number):
        self.cluster_number = cluster_number
        
    def get_data(self,objects_to_cluster):
        self.to_cluster = objects_to_cluster
        
    def calculate_distances(self):
        self.all_distances = numpy.zeros((len(self.to_cluster),len(self.to_cluster)))
        for i in range(len(self.to_cluster)):
            for j in range(len(self.to_cluster)):
                self.all_distances[i,j] = self.distance(self.to_cluster[i],self.to_cluster[j])
                if i%100 == 0:
                    print("i is ",i,end='\r')
        
      
    
    ##function to determine distance between clusters for single linkage returns the indices of the current dendogram level to be lumped together
    ## *****NOTE ***** this method is super slow. I should instead calculate all this distances once and then look them up according to the rules in the following three functions
    def single_linkage(self,dendogram_level):
        unique_cluster_list = numpy.unique(self.dendogram[dendogram_level,:]) #get the list of unique clusters at this level of the dendogram
        #print(unique_cluster_list)
        num_unique_clusters = len(unique_cluster_list) #get the number of unique clusters at this level of the dendogram
        #print(num_unique_clusters)
        temp_distances = numpy.full((num_unique_clusters,num_unique_clusters),1000000000000000) # square matrix to find distance between all existing clusters

        
        for i in range(num_unique_clusters):
            for j in range(num_unique_clusters):
                num_i = len(numpy.where(self.dendogram[dendogram_level,:] == unique_cluster_list[i])[0])
                num_j = len(numpy.where(self.dendogram[dendogram_level,:] == unique_cluster_list[j])[0])

                #print(num_i)
                #print(num_j)
                
                
                if num_i == 1 and num_j == 1:
                    #print(int(unique_cluster_list[i]))
                    #print(numpy.where(self.dendogram[dendogram_level,:]==unique_cluster_list[i])[0][0])
                    #print(self.to_cluster[numpy.where(self.dendogram[dendogram_level,:]==int(unique_cluster_list[i]))[0]])
                    temp_distances[i,j] = self.distance_lookup(self.to_cluster[numpy.where(self.dendogram[dendogram_level,:]==unique_cluster_list[i])[0][0]],self.to_cluster[numpy.where(self.dendogram[dendogram_level,:]==unique_cluster_list[j])[0][0]])
                
                else:
                    short_temp_distances = numpy.zeros((num_i,num_j))
                    

                    if num_i > 1 and num_j == 1:
                        for z in range(num_i):
                                #print(numpy.where(self.dendogram[dendogram_level,:]==unique_cluster_list[i])[0][z])
                                short_temp_distances[z,0] = self.distance_lookup(self.to_cluster[numpy.where(self.dendogram[dendogram_level,:]==unique_cluster_list[i])[0][z]],self.to_cluster[numpy.where(self.dendogram[dendogram_level,:]==unique_cluster_list[j])[0][0]])

                    else:
                        
                        if num_i == 1 and num_j > 1:
                                for zed in range(num_j):
                                    #print(numpy.where(self.dendogram[dendogram_level,:]==unique_cluster_list[j])[0][zed])
                                    short_temp_distances[0,zed] = self.distance_lookup(self.to_cluster[numpy.where(self.dendogram[dendogram_level,:]==unique_cluster_list[i])[0][0]],self.to_cluster[numpy.where(self.dendogram[dendogram_level,:]==unique_cluster_list[j])[0][zed]])

                        else:
                            
                            if num_i > 1 and num_j > 1:
                                for z in range(num_i):
                                    for zed in range(num_j):
                                        #print(numpy.where(self.dendogram[dendogram_level,:]==unique_cluster_list[i])[0][z])
                                        short_temp_distances[z,zed] = self.distance_lookup(self.to_cluster[numpy.where(self.dendogram[dendogram_level,:]==unique_cluster_list[i])[0][z]],self.to_cluster[numpy.where(self.dendogram[dendogram_level,:]==unique_cluster_list[j])[0][zed]])

                
                            
                
                    #print(short_temp_distances)
                    temp_distances[i,j] = numpy.amin(short_temp_distances)
                
                
                if i == j:
                    temp_distances[i,j] = 1000000000000000

               
        #print(temp_distances)      
        #print(numpy.amin(temp_distances))
        first_cluster_index = numpy.where(temp_distances == numpy.amin(temp_distances))[0][0]
        second_cluster_index = numpy.where(temp_distances == numpy.amin(temp_distances))[1][0]
        
        
        return numpy.hstack((numpy.where(self.dendogram[dendogram_level,:] == unique_cluster_list[first_cluster_index])[0],numpy.where(self.dendogram[dendogram_level,:] == unique_cluster_list[second_cluster_index])[0]))       
                
        
        
    #function to determine distance between clusters for complete linkage
    def complete_linkage():
        pass
        
    #function to determine distance between clusters for average linkage
    def average_linkage():
        pass
    
    ##function to look up the distance between two ligands from existing table
    def distance_lookup(self,object_1,object_2):
        return self.all_distances[self.to_cluster.index(object_1),self.to_cluster.index(object_2)]
    
    
    #function to determine distance between two ligands, uses euclidian distance, accepts ligand objects
    def distance(self,object_1,object_2):
        total_distance = len(object_1.OnBits) + len(object_2.OnBits) -  2*numpy.sum(object_1.OnBits == object_2.OnBits)
        euclid_dist = numpy.sqrt(total_distance)
        return euclid_dist
        
    ##function to determine clusters
    def cluster(self):
        self.calculate_distances()
        self.dendogram = numpy.zeros((len(self.to_cluster),len(self.to_cluster)))
        self.dendogram[0,:] = numpy.arange(len(self.to_cluster)) #every object gets its own cluster number
        
        ## goes through 
        for zed in range(1,len(self.to_cluster)-self.cluster_number+1):
            indices_to_lump = self.single_linkage(zed-1) #return an array of indices in the dendogram to lump together
            
            lowest_previous_cluster = numpy.amin(self.dendogram[zed-1,indices_to_lump])
            self.dendogram[zed,:] = self.dendogram[zed-1,:]
            self.dendogram[zed,indices_to_lump] = lowest_previous_cluster

            
            now = datetime.now()
            current_time = now.strftime("%H:%M:%S")
            print("Current Time 2 =", current_time, "zed is ", zed,end='\r')

        
        
    
    
class PartitionClustering():
    ##this class will implement traditional K-means clustering
    def __init__(self,centroid_paradigm,cluster_number):
        self.cluster_number = cluster_number

    def get_data(self,objects_to_cluster):
        self.to_cluster = objects_to_cluster
        self.cluster_assignments = numpy.zeros((len(self.to_cluster),cluster_number))
        self.data_expanded = numpy.zeros((len(self.to_cluster),1024))
        i = 0
        for x in self.to_cluster:
            self.data_expanded[i,x.OnBits] = 1
            i = i + 1
       
        
        
    
    #function to determine distance between two ligands, uses euclidian distance, accepts ligand objects
    def distance(self,object_1,object_2):
        total_distance = len(object_1.OnBits) + len(object_2.OnBits) -  2*numpy.sum(object_1.OnBits == object_2.OnBits)
        euclid_dist = numpy.sqrt(total_distance)
        return euclid_dist
    
    def initialize_centroids(self):
        self.centroids = numpy.zeros((self.cluster_number,len()))
    
    def update_centroids():
        pass
    
    def assign_cluster():
        pass
    
    def is_cluster_over():
        pass
    
    def cluster():
        pass
    
    
    
    

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
