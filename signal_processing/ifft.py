import cmath
from signal_processing import fft

def ifft_2(x):
    """
    Fungsi IFFT untuk menghitung Radix-2.
    Jumlah data berupa pangkat dua.
    """
    N = len(x)
    if N == 1:
        return x
    
    odd_res = ifft_2(x[1::2])
    even_res = ifft_2(x[::2])
    res = [0] * N
    for k in range(N//2):
        W = cmath.exp(1j * 2 * cmath.pi * k / N)
        res[k] = even_res[k] + W*odd_res[k]
        res[k + N//2] = even_res[k] - W*odd_res[k]
    return res

def ifft_3(x):
    """
    Fungsi IFFT untuk menghitung Radix-3.
    Jumlah data berupa pangkat dua.
    """
    N = len(x)
    if N == 1:
        return x
    A_res = ifft_3(x[0::3])
    B_res = ifft_3(x[1::3])
    C_res = ifft_3(x[2::3])
    res = [0] * N
    omega = cmath.exp(1j * 2 * cmath.pi / 3)
    for k in range(N//3):
        W1 = cmath.exp(1j * 2 * cmath.pi * k / N)
        W2 = cmath.exp(1j * 2 * cmath.pi * 2 *  k / N)
        res[k] = A_res[k] + W1 * B_res[k]  + W2 * C_res[k]
        res[k + N//3] = A_res[k] + omega * W1 * B_res[k] + omega**2 * W2 * C_res[k]
        res[k + 2*N//3] = A_res[k] + omega**2 * W1 * B_res[k] + omega * W2 * C_res[k]
    return res
        

def ifft_4(x):
    """
    Fungsi IFFT untuk menghitung Radix-4.
    Jumlah data berupa pangkat dua.
    """
    N = len(x)
    if N == 1:
        return x
    A_res = ifft_4(x[0::4])
    B_res = ifft_4(x[1::4])
    C_res = ifft_4(x[2::4])
    D_res = ifft_4(x[3::4])
    res = [0] * N
    for k in range(N//4):
        W1 = cmath.exp(1j * 2 * cmath.pi * k / N)
        W2 = cmath.exp(1j * 2 * cmath.pi * 2 * k / N)
        W3 = cmath.exp(1j * 2 * cmath.pi * 3 * k / N)
        res[k] = A_res[k] + W1 * B_res[k] + W2 * C_res[k] + W3 * D_res[k]
        res[k+N//4] = A_res[k] + 1j * W1 * B_res[k] - W2 * C_res[k] - 1j * W3 * D_res[k]
        res[k+N//2] = A_res[k] - W1 * B_res[k] + W2 * C_res[k] - W3 * D_res[k]
        res[k+3*N//4] = A_res[k] - 1j * W1 * B_res[k] - W2 * C_res[k] + 1j * W3 * D_res[k]
    return res


def ifft(x, redix = 2):
    """
    Fungsi IFFT untuk menghitung
    """
    if redix == 2:
        output = ifft_2(x)
        return [v/len(x) for v in output]
    elif redix == 3:
        output = ifft_3(x)
        return [v/len(x) for v in output]
    elif redix == 4:
        output = ifft_4(x)
        return [v/len(x) for v in output]
    else:
        raise IndexError("Index Out")


def ifft2d(x):
    """
    
    """
    res = []
    for row in x:
        res.append(ifft(row))

    N = len(res)
    M = len(res[0])

    result = [[0] * M for _ in range(N)]

    for col in range(M):
        column = [res[row][col] for row in range(N)]
        ifft_col = ifft(column)
        for row in range(N):
            result[row][col] = ifft_col[row]

    return result
