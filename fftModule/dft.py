import cmath


def widlet(N):
    """
    Function to making widle matriks
    """
    res = []
    for k in range(N):
        row = []
        for n in range(N):
            W = cmath.exp(-1j * 2 * cmath.pi * k * n / N)
            row.append(W)
        res.append(row)
    return res


def dct(x):
    pass


def dft(x):
    """
    Function to calculation Discrete Fourier Transform
    """
    widle_matriks = widlet(len(x))
    res = [0] * len(x)
    for k in range(len(x)):
        for n in range(len(widle_matriks[0])):
            res[k] += widle_matriks[k][n] * x[k] 
    return res


def dft2d(x, y):
    pass

