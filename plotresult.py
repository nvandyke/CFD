from readgri import readgri
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.collections import PatchCollection
import argparse

#read the state file and return the values
def readresult(filename):
    with open(filename,'r') as f:
        lines = f.readlines()
    arrayholder = []
    for line in lines:
        #we are done when there is an empty line
        line = line.strip()
        if not line:
            continue
        arrayholder.append([float(s) for s in line.split()])
    result = np.array(arrayholder,dtype=float)
    
    return result

#read the error file and return the values
def readerror(filename):
    with open(filename,'r') as f:
        lines = f.readlines()
    vals = []
    for line in lines:
        #we are done on an empty line or non-number
        line = line.strip()
        if line == '':
            break
        line = float(line)
        if not line:
            break
        vals.append(line)
    return vals

#plot vals over the mesh
def plotPDE(mesh,vals):
    
    #set up the plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    patches = []

    #plot mesh with colors
    tris = tri.Triangulation(mesh['V'][:,0],mesh['V'][:,1],mesh['E'])
    tcf = ax.tripcolor(tris,vals,cmap='jet')
    
    #format plot
    ax.set_xlim([np.min(mesh['V'][:,0]),np.max(mesh['V'][:,0])])
    ax.set_ylim([np.min(mesh['V'][:,1]),np.max(mesh['V'][:,1])])
    ax.set_aspect('equal', adjustable='box')
    fig.colorbar(tcf)

#plot the error history
def ploterror(error):

    #create plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    plt.semilogy(error)
    plt.plot([0,len(error)],[1e-7]*2)
    plt.xlabel('Iteration')
    plt.ylabel('Residual')

def main(args):
    mesh = readgri(args.meshfile)
    result = readresult(args.resultfile)
    error = readerror(args.errorfile)

    speed = np.sqrt((result[:,1]/result[:,0])**2 + (result[:,2]/result[:,0])**2)
    y = 1.4
    Pressure = (y-1)*(result[:,3] - .5 * result[:,0]*speed**2)
    c = np.sqrt(y*Pressure/result[:,0])
    Mach = speed/c

    ploterror(error)
    plotPDE(mesh,Mach)
    plt.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        prog='plotResult',
        description='plots solution to gri file and error history')
                                     
    parser.add_argument('meshfile')
    parser.add_argument('resultfile')
    parser.add_argument('errorfile')

    args = parser.parse_args()

    main(args)




