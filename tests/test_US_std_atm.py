import numpy as np
import matplotlib.pyplot as plt

with open('test_US_atm/data.txt', 'r') as f:
    txt = f.read()

txt = txt.split('\n')
headers = txt[0].split(',')
data = np.array([[float(val) for val in line.split(',')] for line in txt[1:] if line != ''])

fig, axes = plt.subplots(nrows=1, ncols=len(headers)-1, sharey=True)

for i, ax in enumerate(axes):
    ax.plot(data[:, i+1], data[:, headers.index('z')] / 1000)
    ax.set_xlabel(headers[i+1])
    if i == 0:
        ax.set_ylabel(headers[0])


plt.show()
