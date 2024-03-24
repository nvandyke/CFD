from readgri import readgri
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.collections import PatchCollection

def readresult(filename):
    with open(filename,'r') as f:
        lines = f.readlines()
    arrayholder = []
    for line in lines:
        line = line.strip()
        if not line:
            continue
        arrayholder.append([float(s) for s in line.split()])
    result = np.array(arrayholder,dtype=float)
    
    return result

def readerror(filename):
    with open(filename,'r') as f:
        lines = f.readlines()
    vals = []
    for line in lines:
        line = line.strip()
        if line == '':
            break
        line = float(line)
        if not line:
            break
        vals.append(line)
    return vals

def normalize(data):
    return (data-np.min(data))/(np.max(data)-np.min(data)) 

def plotPDE(mesh,vals):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    colorMap = plt.get_cmap('jet')
    
    vals = normalize(vals)


    patches = []
    for index, triangle in enumerate(mesh['E']):
        constructed = []
        for vertex in triangle:
            constructed.append(mesh['V'][vertex])
        constructed = np.array(constructed)
        #print(constructed)
        
        tri = plt.Polygon(constructed,color=colorMap(vals[index]))
        patches.append(tri)
        #plt.gca().add_patch(tri)
    ax.add_collection(PatchCollection(patches, match_original=True))
    ax.set_xlim([-2,2])
    ax.set_ylim([-.5,1.5])

def ploterror(error):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    
    plt.semilogy(error)
    plt.plot([0,len(error)],[1e-7]*2)
    plt.xlabel('Iteration')
    plt.ylabel('Residual')

resultfile = 'u.txt'
errorfile = 'e.txt'
meshfile = 'meshes\\kbump2.gri'

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

