import numpy as np
import fftModule as fM
import unittest

class testingfft(unittest.TestCase):
    N = np.linspace(0, 2 * np.pi, 64)
    A = np.sin(N)
    def testfft(self): 
        B = fM.fft.fft(self.A)
        res = np.fft.fft(self.A)
        np.testing.assert_allclose(res, B)

    def testifft(self):
        B = fM.ifft.ifft(self.A)
        res = np.fft.ifft(self.A)
        np.testing.assert_allclose(res, B)
