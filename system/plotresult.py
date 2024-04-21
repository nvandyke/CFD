from readgri import readgri
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.tri as tri
from matplotlib.collections import PatchCollection

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

resultfile = 'bump4solution\\u2.txt'
errorfile = 'bump4solution\\e2.txt'
meshfile = 'meshes\\bump4.gri'

mesh = readgri(meshfile)
result = readresult(resultfile)
error = readerror(errorfile)

speed = np.sqrt((result[:,1]/result[:,0])**2 + (result[:,2]/result[:,0])**2)
y = 1.4
Pressure = (y-1)*(result[:,3] - .5 * result[:,0]*speed**2)
c = np.sqrt(y*Pressure/result[:,0])
Mach = speed/c



#print(mesh)
#print(result)
ploterror(error)
plotPDE(mesh,Mach)
plt.show()

