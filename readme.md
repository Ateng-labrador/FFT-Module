# ProjectFourier - Fast Fourier Transform & Image Compression

ProjectFourier is implementation of various Fast Fourier Transform (FFT)algorithms for signal analysis and 
image compression.This project provides any implementation from any FFT like Cooley-Tukey 
(Radix-2, Radix-3, Radix-4) and Bluestein, including their inverse (IFFT) for data reconstruction.

## Main Feature

- **various**: Radix-2, Radix-3, Radix-4, Cooley-Tukey, and Bluestein
- **DFT & IFFT**: Discrete Fourier transform and invers
- **FFT 2D**: support to transform image (image compression)

## Example

```python
from signal_processing import fft
import numpy as np

# Buat sinyal sample
x = np.array([1, 2, 3, 4, 5, 6, 7, 8], dtype=complex)

# Hitung FFT menggunakan Cooley-Tukey (Radix-2)
X = fft.fft(x)
print("FFT Result:", X)

# Atau gunakan Bluestein (support ukuran arbitrary)
X_bluestein = fft.bluestein(x)
print("Bluestein Result:", X_bluestein)
```

## 📊 Algoritma yang Diimplementasikan

| Algoritma | Kompleksitas | Catatan |
|-----------|--------------|--------|
| DFT | O(n²) | Direct method, slow tapi accurate |
| FFT (Radix-2) | O(n log n) | Cooley-Tukey, untuk n = 2^k |
| FFT (Radix-3/4) | O(n log n) | Cooley-Tukey variant |
| Bluestein | O(n log n) | Arbitrary size, menggunakan Cooley-Tukey internally |
| IFFT | O(n log n) | Inverse transform |
