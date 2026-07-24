import numpy as np
import fftModule as fM
import unittest

class testingbluestein(unittest.TestCase):
    N = np.linspace(0, 2 * np.pi, 1000)
    A = np.sin(N)
    def testbluestein(self):
        B = fM.fft.bluestein(self.A)
        res = np.fft.fft(self.A)
        np.testing.assert_allclose(B, res, atol=1e-7, rtol=1e-7)
