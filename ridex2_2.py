import matplotlib.pyplot as plt
import numpy as np
from signal_processing import fft, ifft
from scipy.signal import find_peaks

# redix 2_2
A1 = 1
A2 = 0.5
A3 = 0.25

f = 100
fs = 1024
N = 1024

t = np.arange(N) / fs

x = (A1 * np.sin(2 * np.pi * f * t) + 
     A2 * np.sin(2 * np.pi * 2 * f * t) + 
     A3 * np.sin(2 * np.pi * 3 * f * t))
freq = np.arange(N) * fs / N


fftredix2 = fft.fft(x)
ifftredix = ifft.ifft(fftredix2)

magnitude = np.abs(fftredix2)
peaks, _ = find_peaks(magnitude, height = 0.1)

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
