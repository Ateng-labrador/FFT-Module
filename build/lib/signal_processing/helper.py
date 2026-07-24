import math


def is_prime(x):
    if x <= 1:
        return False
    for i in range(2, int(math.isqrt(x)) + 1):
        if x % i == 0:
            return False
    return True

def PrimitifRoot(N):
    for g in range(2, N):
        if_primitive = True
        for i in range(1, N - 1):
            if pow(g, i,N) == 1:
                if_primitive = False
                break
        if if_primitive:
            return g
    return None
            

