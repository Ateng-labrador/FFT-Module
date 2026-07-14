from signal_processing import fft, ifft
from PIL import Image
import numpy as np

img = Image.open('sampel.jpeg').convert("L")
img = img.resize((256, 256))
matrix_rgb = np.array(img, dtype=float)

print("=== INPUT ===")
print(f"Min: {matrix_rgb.min()}, Max: {matrix_rgb.max()}, Mean: {matrix_rgb.mean()}")

# Forward FFT
print("\n=== FORWARD FFT ===")
F = fft.bluestein2d(matrix_rgb)
F_array = np.array(F)
F_magnitude = np.abs(F_array)
print(f"Magnitude - Min: {F_magnitude.min()}, Max: {F_magnitude.max()}, Mean: {F_magnitude.mean()}")

# Inverse FFT
print("\n=== INVERSE FFT ===")
iF = ifft.ifft2d(F)
iF_array = np.array([[np.real(iF[i][j]) for j in range(len(iF[0]))] for i in range(len(iF))])
print(f"Min: {iF_array.min()}, Max: {iF_array.max()}, Mean: {iF_array.mean()}")

# Check if values are too small
print(f"\nAda nilai > 1? {(iF_array > 1).sum()}")
print(f"Ada nilai < 0? {(iF_array < 0).sum()}")
print(f"Semua 0? {(iF_array == 0).sum()} dari {iF_array.size}")

# Try scaling
iF_scaled = iF_array / iF_array.max() * 255 if iF_array.max() > 0 else iF_array
print(f"\nSetelah scaling - Min: {iF_scaled.min()}, Max: {iF_scaled.max()}")
