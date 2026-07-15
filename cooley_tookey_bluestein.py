from signal_processing import fft
import matplotlib.pyplot as plt
import time
import numpy as np

samples = 500
result = np.zeros((2, samples))
for i in range(1, samples):
    x_sample = np.random.uniform(-1 * i, i, i) + 1j * np.random.uniform(-1 * i, i, i)
    start_time = time.time()
    fft.fft(x_sample)
    result[0, i] = time.time() - start_time
    start_time = time.time()
    fft.bluestein(x_sample)
    result[1, i] = time.time() - start_time


x_linspace = np.linspace(1, samples, samples)
plt.plot(x_linspace, result[0, :], 'r-', label="cooley-tookey")
plt.plot(x_linspace, result[0, :], 'b-', label="bluestein")
plt.xlabel('n')
plt.ylabel('time')
plt.legend()
plt.show()
