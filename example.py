from signal_processing import fft, ifft
from PIL import Image
import matplotlib.pyplot as plt
import numpy as np


img = Image.open('sampel.jpeg').convert("L")
matrix_rgb = np.array(img, dtype=float)
F = fft.fft2(matrix_rgb)
iF = ifft.ifft2d(F)
image_back = np.real(iF)
image_back = np.clip(image_back, 0, 255)
image_back = image_back.astype(np.uint8)
plt.imshow(image_back, cmap='gray')
plt.show()
