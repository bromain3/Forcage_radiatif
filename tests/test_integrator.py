import numpy as np
import matplotlib.pyplot as plt

with open('test_integrator/data.txt', 'r') as f:
    txt = f.readlines()


data = np.array([[float(val) for val in line.split(',')] for line in txt if line != ''])

plt.plot(data[0,:], data[1,:], label="initial irradiance")
plt.plot(data[0,:], data[2,:], label="irradiance after 1 step of propagation")
plt.legend()
plt.show()

plt.plot(data[0,:], data[1,:], label="initial irradiance")
plt.plot(data[0,:], data[3,:], label="irradiance after full propagation")
plt.legend()
plt.show()

