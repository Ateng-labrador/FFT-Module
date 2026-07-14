import cmath
import math
from signal_processing import ifft


def fft_2(x):
    """
    Fungsi FFT untuk menghitung Radix-2.
    Jumlah data berupa pangkat dua.
    """
    N = len(x)
    if N == 1:
        return x
    odd_res = fft_2(x[1::2])
    even_res = fft_2(x[::2])
    res = [0] * N
    for k in range(N//2):
        W = cmath.exp(-1j * 2 * cmath.pi * k / N)
        res[k] = even_res[k] + W*odd_res[k]
        res[k + N//2] = even_res[k] - W*odd_res[k]
    return res


def fft_3(x):
    """
    Fungsi FFT untuk menghitung Radix-3.
    Jumlah data berupa pangkat dua.
    """
    N = len(x)
    if N == 1:
        return x
    A_res = fft_3(x[0::3])
    B_res = fft_3(x[1::3])
    C_res = fft_3(x[2::3])
    res = [0] * N
    omega = cmath.exp(-1j * 2 * cmath.pi / 3)
    for k in range(N//3):
        W1 = cmath.exp(-1j * 2 * cmath.pi *  k / N)
        W2 = cmath.exp(-1j * 2 * cmath.pi * 2 * k / N)
        res[k] = A_res[k] + W1 * B_res[k] + W2 * C_res[k]
        res[k + N//3] = A_res[k] + omega * W1 * B_res[k] + omega**2 * W2 * C_res[k]
        res[k + 2*N//3] = A_res[k] + omega**2 * W1 * B_res[k] + omega * W2 * C_res[k]
    return res


def fft_4(x):
    """
    Fungsi FFT untuk menghitung Radix-4.
    Jumlah data berupa pangkat dua.
    """
    N = len(x)
    if N == 1:
        return x
    A_res = fft_4(x[0::4])
    B_res = fft_4(x[1::4])
    C_res = fft_4(x[2::4])
    D_res = fft_4(x[3::4])
    res = [0] * N
    for k in range(N//4):
        W1 = cmath.exp(-1j * 2 * cmath.pi * k / N) 
        W2 = cmath.exp(-1j * 2 * cmath.pi * 2 * k / N)
        W3 = cmath.exp(-1j * 2 * cmath.pi * 3 * k / N)
        res[k] = A_res[k] + W1 * B_res[k] + W2 * C_res[k] + W3 * D_res[k]
        res[k+N//4] = A_res[k] - 1j * W1 * B_res[k] - W2 * C_res[k] + 1j * W3 * D_res[k]
        res[k+N//2] = A_res[k] - W1 * B_res[k] + W2 * C_res[k] - W3 * D_res[k]
        res[k+3*N//4] = A_res[k] + 1j * W1 * B_res[k] - W2 * C_res[k] - 1j * W3 * D_res[k]
    return res


def fft(x, radix = 2):
    """
    Fungsi Untuk melakukan kalkulasi fft
    """
    if radix == 2:
        return fft_2(x)
    elif radix == 3:
        return fft_3(x)
    elif radix == 4:
        return fft_4(x)
    else:
        raise IndexError("Index out")


def bluestein(a):
  n = len(a)
  xi = cmath.exp(1j * 2 * cmath.pi / n)

  u = [a[i] * xi**(-(i*i)/2) for i in range(n)]
  v = [xi**((i*i)/2) for i in range(n)]
  v_2 = [xi**(-(i * i)/2) for i in range(n)]

  l = 2 ** math.ceil(math.log2(2 * n - 1))

  u_l = u + [0] * (l - n)

  # Konstruksi v_l yang benar: v[0], v[1], ..., v[n-1], 0s, kemudian v[n-1] conjugate reverse
  aux = v[1:]
  aux.reverse()
  v_l = v + [0] * (l - 2*n + 1) + aux

  dft_u_l = fft(u_l)
  dft_v_l = fft(v_l)

  uv_l = [i * j for i, j in zip(dft_u_l, dft_v_l)]

  # Normalisasi sudah ada di ifft()
  ift_l = ifft.ifft(uv_l, redix=2)
  ift = ift_l[:n]  # Ambil hanya n elements pertama
  dft = [i * j for i, j in zip(v_2, ift)]
  return dft


def bluestein2d(image):
    """
    Fungsi untuk menghitung FFT 2D menggunakan algoritma Bluestein.
    Menerapkan Bluestein per baris dan per kolom.
    """
    res = []
    for row in image:
        if isinstance(row, (list, tuple)):
            res.append(bluestein(list(row)))
        else:
            res.append(bluestein(row.tolist()))
    
    N = len(res)
    M = len(res[0])

    result = [[0+0j] * M for _ in range(N)]

    for col in range(M):
        column = []
        for row in range(N):
            column.append(res[row][col])
        fft_col = bluestein(column)
        for row in range(N):
            result[row][col] = fft_col[row]
    return result


def fft2(image):
    """
    Fungsi untuk menghitung FFT 2D menggunakan algoritma FFT Radix
    """
    res = []
    for row in image:
        res.append(fft(row))
    
    N = len(res)
    M = len(res[0])

    result = [[0] * M for _ in range(N)]

    for col in range(M):
        column = []
        for row in range(N):
            column.append(res[row][col])
        fft_col = fft(column)
        for row in range(N):
            result[row][col] = fft_col[row]
    return result
