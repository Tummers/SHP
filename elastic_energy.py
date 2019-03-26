"""
Code to find the elastic energy of a fluid from its Qtensor file
Qtensor format:
i j k Qxx Qxy Qxz Qyy Qyz vx vy vz qx qy qz q Stressyy Stressyz Stressxz
"""
import numpy as np
import sys

class Qtensor:

    def __init__(self, filename):
        

        self.full_file = np.loadtxt(filename) # reads in file as 2d array

        self.lengthi = int(self.full_file[-1, 0] + 1) # size of system, 0, 0, 31 for 1D
        self.lengthj = int(self.full_file[-1, 1] + 1) # 0, 31, 31 for 2d
        self.lengthk = int(self.full_file[-1, 2] + 1)

        self.init_q_lattices()

        self.err_per_step = 0.00001 # half accuracy in Qtensor file


    def init_q_lattices(self):
        """
        creates 3D q lattices [lengthi, lengthj, lengthk]
        """
        self.Qxx_lattice = np.zeros([self.lengthi, self.lengthj, self.lengthk])
        for n in range(len(self.full_file) - 1):
            i = self.full_file[n, 0]
            j = self.full_file[n, 1]
            k = self.full_file[n, 2]
            qxx_val = self.full_file[n ,3]
                
            self.Qxx_lattice[int(i), int(j), int(k)] = qxx_val

        self.Qxy_lattice = np.zeros([self.lengthi, self.lengthj, self.lengthk])           
        for n in range(len(self.full_file) - 1):
            i = self.full_file[n, 0]
            j = self.full_file[n, 1]
            k = self.full_file[n, 2]
            qxy_val = self.full_file[n, 4]
            
            self.Qxy_lattice[int(i), int(j), int(k)] = qxy_val

        self.Qxz_lattice = np.zeros([self.lengthi, self.lengthj, self.lengthk])
        for n in range(len(self.full_file) - 1):
            i = self.full_file[n, 0]
            j = self.full_file[n, 1]
            k = self.full_file[n, 2]
            qxz_val = self.full_file[n, 5]
            
            self.Qxz_lattice[int(i), int(j), int(k)] = qxz_val

        self.Qyy_lattice = np.zeros([self.lengthi, self.lengthj, self.lengthk])
        for n in range(len(self.full_file) - 1):
            i = self.full_file[n, 0]
            j = self.full_file[n, 1]
            k = self.full_file[n, 2]
            qyy_val = self.full_file[n, 6]
            
            self.Qyy_lattice[int(i), int(j), int(k)] = qyy_val

        self.Qyz_lattice = np.zeros([self.lengthi, self.lengthj, self.lengthk])
        for n in range(len(self.full_file) - 1):
            i = self.full_file[n, 0]
            j = self.full_file[n, 1]
            k = self.full_file[n, 2]
            qyz_val = self.full_file[n, 7]
            
            self.Qyz_lattice[int(i), int(j), int(k)] = qyz_val

        self.Qzz_lattice = -self.Qxx_lattice - self.Qyy_lattice # **subtraction or multiplication?**


#-----------------------------------------------------------------------------------------------
    # summation dz components
    def dzQxx(self, i, j, k):
        numerator = self.Qxx_lattice[i, j, k+1] - self.Qxx_lattice[i, j, k-1]
        return numerator / 2

    def dzQxy(self, i, j, k):
        numerator = self.Qxy_lattice[i, j, k+1] - self.Qxy_lattice[i, j, k-1]
        return numerator / 2
   
    def dzQxz(self, i, j, k):
        numerator = self.Qxz_lattice[i, j, k+1] - self.Qxz_lattice[i, j, k-1]
        return numerator / 2

    def dzQyy(self, i, j, k):
        numerator = self.Qyy_lattice[i, j, k+1] - self.Qyy_lattice[i, j, k-1]
        return numerator / 2
    
    def dzQyz(self, i, j, k):
        numerator = self.Qyz_lattice[i, j, k+1] - self.Qyz_lattice[i, j, k-1]
        return numerator / 2

    def dzQzz(self, i, j, k):
        numerator = self.Qzz_lattice[i, j, k+1] - self.Qzz_lattice[i, j, k-1]
        return numerator / 2


    # dy components
    def dyQxx(self, i, j, k):
        numerator = self.Qxx_lattice[i, j+1, k] - self.Qxx_lattice[i, j-1, k]
        return numerator / 2
        
    def dyQxy(self, i, j, k):
        numerator = self.Qxy_lattice[i, j+1, k] - self.Qxy_lattice[i, j-1, k]
        return numerator / 2
            
    def dyQxz(self, i, j, k):
        numerator = self.Qxz_lattice[i, j+1, k] - self.Qxz_lattice[i, j-1, k]
        return numerator / 2
        
    def dyQyy(self, i, j, k):
        numerator = self.Qyy_lattice[i, j+1, k] - self.Qyy_lattice[i, j-1, k]
        return numerator / 2

    def dyQyz(self, i, j, k):
        numerator = self.Qyz_lattice[i, j+1, k] - self.Qyz_lattice[i, j-1, k]
        return numerator / 2
             
    def dyQzz(self, i, j, k):
        numerator = self.Qzz_lattice[i, j+1, k] - self.Qzz_lattice[i, j-1, k]
        return numerator / 2

#---------------------------------------------------------------------------------------------------
    def summation_1d(self):
        """
        summation for 1 dimension
        """
        tot = 0
        
        i = 0
        j = 0
        err = 0
        for k in range(1, self.lengthk - 1): # range excludes edge effects
            term1 = self.dzQxx(i, j, k) ** 2
            term2 = 2 * (self.dzQxy(i, j, k) ** 2)
            term3 = 2 * (self.dzQxz(i, j, k) ** 2)
            term4 = self.dzQyy(i, j, k) ** 2
            term5 = 2 * (self.dzQyz(i, j, k) ** 2)
            term6 = self.dzQzz(i, j, k) ** 2 
            
            temporary_tot = term1 + term2 + term3 + term4 + term5 + term6
            err = np.sqrt((err ** 2) + (self.err_per_step ** 2))
            tot += temporary_tot
        
        infile = np.loadtxt("in.txt")
        activity = infile[0]
        
        self.add_to_file(np.array([activity, tot, err]))
        
        return tot, err


    def summation_2d(self):
        """
        summation for 2d
        """
        tot = 0
        
        i = 0
        err = 0
        for j in range(1, self.lengthj - 1):
            for k in range(1, self.lengthk - 1):
                term1 = self.dzQxx(i, j, k) ** 2
                term2 = 2 * (self.dzQxy(i, j, k) ** 2)
                term3 = 2 * (self.dzQxz(i, j, k) ** 2)
                term4 = self.dzQyy(i, j, k) ** 2
                term5 = 2 * (self.dzQyz(i, j, k) ** 2)
                term6 = self.dzQzz(i, j, k) ** 2 
                
                term7 = self.dyQxx(i, j, k) ** 2
                term8 = 2 *(self.dyQxy(i, j, k) ** 2)
                term9 = 2 *( self.dyQxz(i, j, k) ** 2)
                term10 = self.dyQyy(i, j, k) ** 2
                term11 = 2 * (self.dyQyz(i, j, k) ** 2)
                term12 = self.dyQzz(i, j, k) ** 2
                
                temp_tot1 = term1 + term2 + term3 + term4 + term5 + term6
                temp_tot2 = term7 + term8 + term9 + term10 + term11 + term12
                temp_tot = temp_tot1 + temp_tot2

                err = np.sqrt((err ** 2) + (self.err_per_step ** 2))

                tot += temp_tot

        infile = np.loadtxt("in.txt")
        activity = infile[0]
        
        self.add_to_file(np.array([activity, tot, err]))
        
        return tot, err


    def add_to_file(self, line):
        """
        writes a line to a file
        """
        f = open("elastic_energies.txt", "a")
        np.savetxt(f, line, newline=" ")
        f.write("\n")
        f.close()


def input_increment(name, increment):
    """
    increments a value in the in.txt file
    """
    if(name == "activity"):
        index = 0
    elif(name == "bodyforce"):
        index = 1
    elif(name == "shearv"):
        index = 2
    else:
        print("invalid variable name")

    array = np.loadtxt("in.txt")
    array[index] += increment
    print("Current activity is: %.5f" %(array[index]))
    f = open("in.txt", "w")
    for i in range(3):
        f.write(str(array[i]) + " ")


def main():
    
    filename = sys.argv[1]
    dimensions = int(sys.argv[2])
    increment_name = sys.argv[3]
    increment_value = float(sys.argv[4])
    qdata = Qtensor(filename)

    if(dimensions == 1): 
        elastic_energy, err = qdata.summation_1d()
    elif(dimensions == 2): 
        elastic_energy, err = qdata.summation_2d()
    else:
        print("Input as: 'python3 elastic_energy.py <filename> <number of dimensions> <increment_name> <increment_value>'")
        exit()
    
        
    input_increment(increment_name, increment_value)
    
    # print("Elastic energy: %f, Error: %f" %(elastic_energy, err))

main()
