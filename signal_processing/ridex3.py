import matplotlib.pyplot as plt
import numpy as np
import fft
import ifft

# redix 3
t = np.linspace(0, 2*np.pi, 81)
x = np.random.rand(len(t))
fftredix3 = fft.fft(x, radix=3)
ifftredix3 = ifft.ifft(fftredix3, redix=3)

fig, axis = plt.subplots(2, 2, figsize=(10, 8))
axis[0, 0].plot(t, x)
axis[0, 0].set_title("Sinyal Random")


axis[0, 1].plot(np.real(fftredix3))
axis[0, 1].set_title("Hasil FFT")


axis[1, 0].plot(t, np.abs(ifftredix3))
axis[1, 0].set_title("Hasil IFFT")


axis[1, 1].plot(t, x, label="Original")
axis[1, 1].plot(t , np.abs(ifftredix3), "--" ,label="IFFT")
axis[1, 1].set_title("Fitting")
axis[1, 1].legend()

plt.show()

