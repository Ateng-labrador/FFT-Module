import matplotlib.pyplot as plt
import numpy as np
from signal_processing import fft, ifft
from scipy.signal import find_peaks

A = 1
fs = 1024
f = 100
N = 1024

t = np.arange(N) / fs

x = A * np.sin(2 * np.pi * f * t)
freq = np.arange(N) * fs / N

fftredix2 = fft.bluestein(x)
ifftredix2 = ifft.ifft(fftredix2)

magnitude = 2 * np.abs(fftredix2) / N
peaks, _ = find_peaks(magnitude, height=0.1)

plt.plot(freq, magnitude)
plt.scatter(freq[peaks], magnitude[peaks], color='red')

for p in peaks:
    plt.annotate(
        f"{freq[p]:.2f}",
        xy=(freq[p], magnitude[p]),
        xytext=(0, 8),
        textcoords="offset points",
        ha="center"
    )

plt.show()

