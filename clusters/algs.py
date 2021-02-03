class HierarchicalClustering():
    import numpy
    def __init__(self,distance_paradigm,cluster_number):
        self.dist_prdgm = distance_paradigm
        self.dendogram  = [] #somethings a matric maybe?
        self.cluster_number = cluster_number
        
    def update_cluster(self,cluster_number):
        self.cluster_number = cluster_number
        
    def get_data(self,objects_to_cluster):
        self.to_cluster = objects_to_cluster
        
    def calculate_distances(self):
        pass
        
      
    
    ##function to determine distance between clusters for single linkage returns the indices of the current dendogram to be lumped together
    ## *****NOTE ***** this method is super slow. I should instead calculate all this distances once and then look them up according to the rules in the following three functions
    def single_linkage(self):
        import numpy
        temp_distances = numpy.zeros((len(self.dendogram),len(self.dendogram))) # square matrix to find distance between all existing clusters
        for i in range(len(self.dendogram)):
            for j in range(len(self.dendogram)):
                num_i = 1
                num_j = 1
                if type(self.dendogram[i]) == list:
                    num_i = (len(self.dendogram[i]))
                             
                if type(self.dendogram[j]) == list:
                    num_j = (len(self.dendogram[j]))

                
                short_temp_distances = numpy.zeros((num_i,num_j))
                
                
                if num_i > 1 and num_j > 1:
                    for z in range(num_i):
                        for zed in range(num_j):
                            short_temp_distances[z,zed] = self.distance(self.dendogram[i][z],self.dendogram[j][zed])
                            
                if num_i > 1 and num_j == 1:
                    for z in range(num_i):
                        for zed in range(num_j):
                            short_temp_distances[z,zed] = self.distance(self.dendogram[i][z],self.dendogram[j])   
                
                if num_i == 1 and num_j > 1:
                    for z in range(num_i):
                        for zed in range(num_j):
                            short_temp_distances[z,zed] = self.distance(self.dendogram[i],self.dendogram[j][zed])
                            
                if num_i == 1 and num_j == 1:
                    for z in range(num_i):
                        for zed in range(num_j):
                            short_temp_distances[z,zed] = self.distance(self.dendogram[i],self.dendogram[j])
                            
                if i == j:
                    temp_distances[i,j] = 1000000000000000
                    
                if i != j:
                    temp_distances[i,j] = numpy.amin(short_temp_distances)
                
                print("i is ",i,end='\r')
               
                
        pair_first_index = numpy.where(temp_distances == numpy.amin(temp_distances))[0][0]
        pair_second_index = numpy.where(temp_distances == numpy.amin(temp_distances))[1][0]
        return [pair_first_index, pair_second_index]       
                
        
        pass
    #function to determine distance between clusters for complete linkage
    def complete_linkage():
        pass
        
    #function to determine distance between clusters for average linkage
    def average_linkage():
        pass
    
    
    #function to determine distance between two ligands, uses euclidian distance, accepts ligand objects
    def distance(self,object_1,object_2):
        import numpy
        total_distance = len(object_1.OnBits) + len(object_2.OnBits) -  2*numpy.sum(object_1.OnBits == object_2.OnBits)
        euclid_dist = numpy.sqrt(total_distance)
        return euclid_dist
        
    ##function to determine clusters
    def cluster(self):
        
        ##the first round is every object gets its own cluster
        for z in range(len(self.to_cluster)):
            self.dendogram.append(self.to_cluster[z])
            
        
        ## goes through 
        for zed in range(1,len(self.to_cluster)-self.cluster_number):
            i = self.single_linkage()[0]
            j = self.single_linkage()[1]
            
            #combines objects to be clustered in this interation in a temporary list
            if type(self.dendogram[i]) == list and type(self.dendogram[j]) == list:
                dendogram_temp = self.dendogram[i] + self.dendogram[j]
            if type(self.dendogram[i]) != list and type(self.dendogram[j]) == list:
                dendogram_temp = [self.dendogram[i]] + self.dendogram[j]
            if type(self.dendogram[i]) == list and type(self.dendogram[j]) != list:
                dendogram_temp = self.dendogram[i] + [self.dendogram[j]]
            if type(self.dendogram[i]) != list and type(self.dendogram[j]) != list:
                dendogram_temp = [self.dendogram[i]] + [self.dendogram[j]]
            
            #first creates a new list with out the new cluster
            my_new_dendogram = []
            for k in range(len(self.dendogram)):
                if j != k and i != k:
                    my_new_dendogram.append(self.dendogram[k])

            my_new_dendogram.append(my_dendogram_temp) #adds the new cluster back
            self.dendogram = my_new_dendogram #updates the dendogram cluster list in the object
            print("zed is ",zed)

        
        
    
    
class PartitionClustering():
    ##this class will implement K-means clustering
    def __init__(self,centroid_paradigm,cluster_number):
        self.hmm = cluster_number

    ##needs a class to init
    
    
    
    

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
