import numpy as np
import matplotlib.pyplot as plt

with open('test_data_handler/tau.txt', 'r') as f:
    txt = f.readlines()


data = np.array([[float(val) for val in line.split(',')] for line in txt if line != ''])

plt.plot(data[:,0], data[:,1])
plt.yscale('log')
plt.show()

plt.plot(data[:,0], data[:,2])
plt.show()

plt.plot(data[:,0], data[:,3])
plt.show()
