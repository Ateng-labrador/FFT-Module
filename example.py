from signal_processing import fft, ifft
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np

# Load dan resize gambar
img = Image.open('sampel.jpeg').convert("L")
img = img.resize((256, 256))
matrix_rgb = np.array(img, dtype=float)

print("Starting FFT 2D Compression...")

# Forward FFT 2D dengan Bluestein
F = fft.bluestein2d(matrix_rgb)
F_array = np.array(F)

# Hitung magnitude untuk compression
F_magnitude = np.abs(F_array)
total_coeffs = F_magnitude.size

# COMPRESSION: Hapus low-magnitude frequencies
# Ubah threshold untuk kontrol compression ratio
compression_ratios = [0.50, 0.75, 0.90]  # 50%, 75%, 90% compression

fig, axes = plt.subplots(2, 2, figsize=(12, 12))

# Original image
axes[0, 0].imshow(matrix_rgb, cmap='gray')
axes[0, 0].set_title('Original Image')
axes[0, 0].axis('off')

# Reconstruct tanpa compression
iF = ifft.ifft2d(F)
image_back = np.array([[np.real(iF[i][j]) for j in range(len(iF[0]))] for i in range(len(iF))])
image_back = np.clip(image_back, 0, 255).astype(np.uint8)

axes[0, 1].imshow(image_back, cmap='gray')
axes[0, 1].set_title('Reconstructed (No Compression)')
axes[0, 1].axis('off')
error_0 = np.mean((matrix_rgb - image_back) ** 2)
axes[0, 1].text(0.5, -0.1, f'RMSE: {error_0:.2f}', ha='center', transform=axes[0, 1].transAxes)

# Compression dengan threshold berbeda
for idx, ratio in enumerate(compression_ratios[:2]):
    # Set threshold untuk mencapai compression ratio
    threshold = np.percentile(F_magnitude, ratio * 100)
    
    # Copy FFT dan hapus low-magnitude values
    F_compressed = F_array.copy()
    mask = F_magnitude < threshold
    F_compressed[mask] = 0
    
    removed = mask.sum()
    print(f"\nCompression {ratio*100:.0f}%:")
    print(f"  Coefficients removed: {removed}/{total_coeffs}")
    print(f"  Threshold: {threshold:.2f}")
    
    # Inverse FFT
    iF_comp = ifft.ifft2d(F_compressed.tolist())
    image_comp = np.array([[np.real(iF_comp[i][j]) for j in range(len(iF_comp[0]))] for i in range(len(iF_comp))])
    image_comp = np.clip(image_comp, 0, 255).astype(np.uint8)
    
    # Calculate RMSE
    rmse = np.sqrt(np.mean((matrix_rgb - image_comp) ** 2))
    
    axes[1, idx].imshow(image_comp, cmap='gray')
    axes[1, idx].set_title(f'Compressed ({ratio*100:.0f}%)')
    axes[1, idx].axis('off')
    axes[1, idx].text(0.5, -0.1, f'RMSE: {rmse:.2f}', ha='center', transform=axes[1, idx].transAxes)

plt.tight_layout()
plt.savefig('compression_results.png', dpi=100, bbox_inches='tight')
print("\nResults saved to: compression_results.png")
plt.show()

