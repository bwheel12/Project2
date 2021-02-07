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
                    print("i is ",i)
        
      
    
    ##function to determine distance between clusters for single linkage returns the indices of the current dendogram to be lumped together
    ## *****NOTE ***** this method is super slow. I should instead calculate all this distances once and then look them up according to the rules in the following three functions
    def single_linkage(self):
        temp_distances = numpy.zeros((len(self.dendogram[self.dendogram_level]),len(self.dendogram[self.dendogram_level]))) # square matrix to find distance between all existing clusters
        for i in range(len(self.dendogram[self.dendogram_level])):
            for j in range(len(self.dendogram[self.dendogram_level])):
                num_i = 1
                num_j = 1
                if type(self.dendogram[self.dendogram_level][i]) == list:
                    num_i = (len(self.dendogram[self.dendogram_level][i]))
                             
                if type(self.dendogram[self.dendogram_level][j]) == list:
                    num_j = (len(self.dendogram[self.dendogram_level][j]))

                
                short_temp_distances = numpy.zeros((num_i,num_j))
                
                if num_i == 1 and num_j == 1:
                            short_temp_distances[0,0] = self.distance_lookup(self.dendogram[self.dendogram_level][i],self.dendogram[self.dendogram_level][j])
                
                else:
                
                    

                    if num_i > 1 and num_j == 1:
                        for z in range(num_i):
                                short_temp_distances[z,0] = self.distance_lookup(self.dendogram[self.dendogram_level][i][z],self.dendogram[self.dendogram_level][j])   

                    else:
                        
                        if num_i == 1 and num_j > 1:
                                for zed in range(num_j):
                                    short_temp_distances[0,zed] = self.distance_lookup(self.dendogram[self.dendogram_level][i],self.dendogram[self.dendogram_level][j][zed])

                        else:
                            
                            if num_i > 1 and num_j > 1:
                                for z in range(num_i):
                                    for zed in range(num_j):
                                        short_temp_distances[z,zed] = self.distance_lookup(self.dendogram[self.dendogram_level][i][z],self.dendogram[self.dendogram_level][j][zed])

                
                            
                
                    
                if i != j:
                    temp_distances[i,j] = numpy.amin(short_temp_distances)
                
                
                else:
                    if i == j:
                        temp_distances[i,j] = 1000000000000000
                
                #if i%50 == 0:
                    #print("Single Linkage is ",i,end='\r')
               
        #print(temp_distances)\
        print(numpy.amin(temp_distances))
        pair_first_index = numpy.where(temp_distances == numpy.amin(temp_distances))[0][0]
        pair_second_index = numpy.where(temp_distances == numpy.amin(temp_distances))[1][0]
        return [pair_first_index, pair_second_index]       
                
        
        
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
        
        ##the first round is every object gets its own cluster
        for z in range(len(self.to_cluster)):
            self.dendogram.append(self.to_cluster[z])
        self.dendogram = [self.dendogram]
        self.dendogram_level = 0
        
        ## goes through 
        for zed in range(1,len(self.to_cluster)-self.cluster_number+1):
            i = self.single_linkage()[0]
            j = self.single_linkage()[1]
            
            
            now = datetime.now()

            current_time = now.strftime("%H:%M:%S")
            print("Current Time 1 =", current_time)
            
            #combines objects to be clustered in this interation in a temporary list
            if type(self.dendogram[self.dendogram_level][i]) == list and type(self.dendogram[self.dendogram_level][j]) == list:
                dendogram_temp = self.dendogram[self.dendogram_level][i] + self.dendogram[self.dendogram_level][j]
            if type(self.dendogram[self.dendogram_level][i]) != list and type(self.dendogram[self.dendogram_level][j]) == list:
                dendogram_temp = [self.dendogram[self.dendogram_level][i]] + self.dendogram[self.dendogram_level][j]
            if type(self.dendogram[self.dendogram_level][i]) == list and type(self.dendogram[self.dendogram_level][j]) != list:
                dendogram_temp = self.dendogram[self.dendogram_level][i] + [self.dendogram[self.dendogram_level][j]]
            if type(self.dendogram[self.dendogram_level][i]) != list and type(self.dendogram[self.dendogram_level][j]) != list:
                dendogram_temp = [self.dendogram[self.dendogram_level][i]] + [self.dendogram[self.dendogram_level][j]]
            
            #first creates a new list with out the new cluster
            my_new_dendogram = []
            for k in range(len(self.dendogram[self.dendogram_level])):
                if j != k and i != k:
                    my_new_dendogram.append(self.dendogram[self.dendogram_level][k])

            my_new_dendogram.append(dendogram_temp) #adds the new cluster back
            self.dendogram.append(my_new_dendogram) #add the new dendogram cluster list ot the last entry in the existing dendogram object 
            self.dendogram_level = self.dendogram_level + 1
            now = datetime.now()

            current_time = now.strftime("%H:%M:%S")
            print("Current Time 2 =", current_time)
            print("zed is ",zed,end='\r')

        
        
    
    
class PartitionClustering():
    ##this class will implement K-means clustering
    def __init__(self,centroid_paradigm,cluster_number):
        self.cluster_number = cluster_number

    def get_data(self,objects_to_cluster):
        self.to_cluster = objects_to_cluster
    
    #function to determine distance between two ligands, uses euclidian distance, accepts ligand objects
    def distance(self,object_1,object_2):
        total_distance = len(object_1.OnBits) + len(object_2.OnBits) -  2*numpy.sum(object_1.OnBits == object_2.OnBits)
        euclid_dist = numpy.sqrt(total_distance)
        return euclid_dist
    
    def initialize_centroids():
        pass
    
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
