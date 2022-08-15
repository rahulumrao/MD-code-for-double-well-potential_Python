#read the data from file
#!/usr/bin/python
'''
import matplotlib.pyplot as plt
import csv
#
X, Y = [], []
with open ('noise.txt', 'r') as csvfile:
    plots = csv.reader(csvfile, delimiter=' ')
    for row in plots:
        X.append(row[0])
        Y.append(row[1])
#
plt.plot(X, Y, Label="noise")
plt.show()
'''
'''
import matplotlib.pyplot as plt
import numpy as np
#
x, y = np.loadtxt('noise.txt', delimiter=' ', unpack=True)
plt.plot(x,y, label='Loaded from file!')

#plt.xlabel('x')
#plt.ylabel('y')
#plt.title('Interesting Graph\nCheck it out')
#plt.legend()
plt.show()
'''
'''
import matplotlib.pyplot as plt
import numpy as np
a = np.loadtxt("ener.txt", unpack=True)
plt.plot(a[0], a[6], label = "data0");
plt.legend();
plt.show();
#print(a[0])
'''
'''
from matplotlib import pyplot as plt;
    from pylab import genfromtxt;  
    mat0 = genfromtxt("Temp.txt");
    #mat1 = genfromtxt("ener.txt");
    plt.plot(mat0[:,0], mat0[:,1], label = "data0");
    #plt.plot(mat1[:,0], mat1[:,1], label = "data1");
    plt.legend();
    plt.show();
    #plot()
'''
ener = open("ener.txt","r")
def plot(file):
    import matplotlib.pyplot as plt
    import numpy as np
    a = np.loadtxt(file, unpack=True)
    plt.plot(a[0], a[6], label = "Langevin Thermostat (NVT)");
    plt.xlabel('MD Steps')
    plt.ylabel('Temperature')
    plt.legend();
    plt.show();
plot('ener.txt')
plot('noise.txt')
