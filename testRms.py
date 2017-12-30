#!/users/rprestag/venv/bin/python

import numpy as np
import math
import pylab as plt

# array sizes
size = 1024*16

# input mean power
inp = 20.0

# deduced rms for each of real and imaginary
rms = math.sqrt(inp/2.0)

# create two arrays for the real and imaginary components

real = np.random.normal(loc=0.0, scale=rms, size=size)
imag = np.random.normal(loc=0.0, scale=rms, size=size)

power = real**2 + imag**2

plt.plot(real, 'r-')
plt.plot(imag, 'b-')
plt.plot(power,'g-')
plt.show()
print power.mean()
