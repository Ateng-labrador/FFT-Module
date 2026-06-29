import matplotlib.pyplot as plt
import numpy as np
import fft
import ifft

# redix 4
t = np.linspace(0, 2*np.pi, 64)
x = np.random.rand(len(t))
fftredix4 = fft.fft(x, radix = 4)
ifftredix4 = ifft.ifft(fftredix4, redix=4)

fig, axis = plt.subplots(2, 2, figsize=(10, 8))
axis[0, 0].plot(t, x)
axis[0, 0].set_title("Sinyal Random")


axis[0, 1].plot(np.real(fftredix4))
axis[0, 1].set_title("Hasil FFT")


axis[1, 0].plot(t, np.abs(ifftredix4))
axis[1, 0].set_title("Hasil IFFT")


axis[1, 1].plot(t, x, label="Original")
axis[1, 1].plot(t , np.abs(ifftredix4), "--" ,label="IFFT")
axis[1, 1].set_title("Fitting")
axis[1, 1].legend()

plt.show()

