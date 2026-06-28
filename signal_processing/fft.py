import cmath

def fft(x):
    N = len(x)
    
    if N == 1:
        return x
    odd_res = fft(x[1::2])
    even_res = fft(x[::2])
    res = [0] * N
    for k in range(N//2):
        W = cmath.exp(-1j * 2 * cmath.pi * k / N)
        res[k] = even_res[k] + W*odd_res[k]
        res[k + N//2] = even_res[k] - W*odd_res[k]
    return res


