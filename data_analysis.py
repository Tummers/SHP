"""
Python programme to analyse data produced by active_rheology.c
"""

import numpy as np
import matplotlib.pyplot as plt
import sys


class Qtensor:

    def __init__(self, data_file, input_file):
        """
        object comprised of array length n with 18 columns

        Input file has columns:
        0 1 2 3   4   5   6   7   8  9  10 11 12 13 14 15       16       17
        i j k Qxx Qxy Qxz Qyy Qyz vx vy vz qx qy qz q  Stressyy Stressyz Stressxz
        """

        self.full_data = np.loadtxt(data_file)
        self.length = len(self.full_data)
        self.input_data = np.loadtxt(input_file)

        self.i = self.full_data[:, 0] # dividing data into its constituent parts
        self.j = self.full_data[:, 1]
        self.k = self.full_data[:, 2]
        self.Qxx = self.full_data[:, 3]
        self.Qxy = self.full_data[:, 4]
        self.Qxz = self.full_data[:, 5]
        self.Qyy = self.full_data[:, 6]
        self.Qyz = self.full_data[:, 7]
        self.vx = self.full_data[:, 8]
        self.vy = self.full_data[:, 9]
        self.vz = self.full_data[:, 10]
        self.qx = self.full_data[:, 11]
        self.qy = self.full_data[:, 12]
        self.qz = self.full_data[:, 13]
        self.q = self.full_data[:, 14]
        self.Stressyy = self.full_data[:, 15]
        self.Stressyz = self.full_data[:, 16]
        self.Stressxz = self.full_data[:, 17]

        self.activity = self.input_data[0]
        self.body_force = self.input_data[1]
        self.shear_v = self.input_data[2]


class directorField:

    def __init__(self, data_file, input_file):
        """
        object comprised of array made from director field file

        Input file has columns:
        0 1 2 3  4  5
        i j k nx ny nz
        """

        self.full_data = np.loadtxt(data_file)
        self.length = len(self.full_data)
        self.input_data = np.loadtxt(input_file)

        self.i = self.full_data[:, 0]
        self.j = self.full_data[:, 1]
        self.k = self.full_data[:, 2]
        self.nx = self.full_data[:, 3]
        self.ny = self.full_data[:, 4]
        self.nz = self.full_data[:, 5]

        self.activity = self.input_data[0]
        self.body_force = self.input_data[1]
        self.shear_v = self.input_data[2]
# -----------------------------------------------------------------------------

def plot(xs, ys, title="", xname="", yname=""):
    """
    Plots xs against ys in matplotlib
    """
    plt.plot(xs, ys)
    plt.title(title)
    plt.xlabel(xname)
    plt.ylabel(yname)
    plt.show()

def add_to_file(line):
    f = open("visc_with_n.txt", "a")
    np.savetxt(f, line, newline=" ")
    f.write("\n")
    f.close()


def activity_viscosity_write_1d(data, field):
    """
    viscosity approaches 0 around 0.005 / 0.006 activity
    """
    temp_array = np.zeros(len(data.Stressyz))
    len_Stressyz = len(data.Stressyz)

    shear_rate = 2 * data.shear_v

    gamma = shear_rate / field.k[field.length - 1]

    avg_Stress = np.average(data.Stressyz)

    pi_yz = avg_Stress

    visc = pi_yz / gamma

    for i in range(len_Stressyz):
        temp_array[i] = (data.Stressyz[i] - avg_Stress) ** 2

    error = np.sqrt(np.sum(temp_array) / (len_Stressyz - 1))

    add_to_file(np.array([data.activity, visc, error]))


def activity_viscosity_write_2d(data, field):
    """
    viscosity approaches 0 around 0.005 / 0.006 activity
    """
    temp_array = np.zeros(len(data.Stressyz))
    len_Stressyz = len(data.Stressyz)

    shear_rate = 2 * data.shear_v

    gamma = shear_rate / field.k[field.length - 1]**2

    avg_Stress = np.average(data.Stressyz)

    pi_yz = avg_Stress / (field.k[field.length - 1])

    visc = pi_yz / gamma

    for i in range(len_Stressyz):
        temp_array[i] = (data.Stressyz[i] - avg_Stress) ** 2

    error = np.sqrt(np.sum(temp_array) / (len_Stressyz - 1))

    add_to_file(np.array([data.activity, visc, error]))


def activity_viscosity_plot():

    array = np.loadtxt("visc_with_n.txt")
    xs = array[:, 0]
    ys = array[:, 1]

    plot(xs, ys, "Viscosity against Activity", "Activity", "Viscosity")


def input_increment(name, increment):
    """
    increments a value of the in.txt file
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

    data = Qtensor("Qtensor.txt", "in.txt")
    field = directorField("directorField.txt", "in.txt")
    mode = sys.argv[1]

    if(mode == "write_1d"):
        name = sys.argv[2]
        increment = float(sys.argv[3])
        activity_viscosity_write_1d(data, field)
        input_increment(name, increment)

    if(mode == "write_2d"):
        name = sys.argv[2]
        increment = float(sys.argv[3])
        activity_viscosity_write_2d(data, field)
        input_increment(name, increment)

    if(mode == "plot"):
        activity_viscosity_plot()

main()
