def egcd(a, b):
    if b == 0:
        return a, 1, 0
    g, y, x = egcd(b, a % b)
    return g, x, y - (a // b) * x

def mod_inverse(a, m):
    g, x, _ = egcd(a, m)
    if g != 1:
        raise ValueError("Оберненого не існує")
    return x % m

def montgomery_setup(N):
    k = N.bit_length()
    R = 1 << (k + 1)  # R > N
    R_inv = mod_inverse(R, N)
    N_prime = (-mod_inverse(N, R)) % R
    return R, R_inv, N_prime

def REDC(T, N, R, N_prime):
    # Montgomery reduction
    m = ((T % R) * N_prime) % R
    t = (T + m * N) // R
    if t >= N:
        t -= N
    return t

def to_montgomery_space(x, N, R):
    return (x * R) % N

def from_montgomery_space(xM, N, R, N_prime):
    return REDC(xM, N, R, N_prime)

def montgomery_multiplication(aM, bM, N, R, N_prime):
    return REDC(aM * bM, N, R, N_prime)

def montgomery_pow_binary(base, exponent, N):
    # Бінарне піднесення до степеня в просторі Монтгомері
    R, R_inv, N_prime = montgomery_setup(N)
    baseM = to_montgomery_space(base, N, R)
    resultM = to_montgomery_space(1, N, R)

    for bit in bin(exponent)[2:]:
        resultM = montgomery_multiplication(resultM, resultM, N, R, N_prime)
        if bit == '1':
            resultM = montgomery_multiplication(resultM, baseM, N, R, N_prime)

    return from_montgomery_space(resultM, N, R, N_prime)

# Приклад
N = 101
base = 7
exp = 13

# Wolfram alpha видав 75
print(f"{base}^{exp} mod {N} = {montgomery_pow_binary(base, exp, N)}")