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


def normalize(data):
    return (data-np.min(data))/(np.max(data)-np.min(data)) 


resultfile= 'u2.txt'
meshfile = 'bump2.gri'
#meshfile = 'test.gri'

mesh = readgri(meshfile)
result = readresult(resultfile)
speed = np.sqrt((result[:,1]/result[:,0])**2 + (result[:,2]/result[:,0])**2)
y = 1.4
Pressure = (y-1)*(result[:,3] - .5 * result[:,0]*speed**2)
c = np.sqrt(y*Pressure/result[:,0])
Mach = speed/c


vals = normalize(Mach)

#print(mesh)
#print(result)
fig = plt.figure()
ax = fig.add_subplot(111)
colorMap = plt.get_cmap('jet')


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
plt.show()
