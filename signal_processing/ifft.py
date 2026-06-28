import cmath

def ifft_(x):
    N = len(x)
    if N == 1:
        return x
    
    odd_res = ifft_(x[1::2])
    even_res = ifft_(x[::2])
    res = [0] * N
    for k in range(N//2):
        W = cmath.exp(1j * 2 * cmath.pi * k / N)
        res[k] = (even_res[k] + W*odd_res[k]) / N
        res[k + N//2] = (even_res[k] - W*odd_res[k]) / N
    return res

def ifft(x):
    output = ifft_(x)
    return [v/len(x) for v in output]
