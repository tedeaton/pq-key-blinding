import random
import hashlib
import string
import time
import csv

## PARAMETERS

p = 35184372067549
n = 512
k = 4
l = 4
d = 7
eta = 7
beta = 644
n1c = 46
gamma1 = 905679
gamma2 = 905679
seedbytes = 32

## RING SETUP

P.<t> = PolynomialRing(Integers(p))
Rq.<x> = QuotientRing(P, P.ideal(t^n + 1))

## POLY/VEC/MAT INFRASTRUCTURE

def apply_poly(f, *args):
    return Rq(sum([f(*(int(coeff) for coeff in coeffs), args[-1])*x^i for i, coeffs in enumerate(zip(*(poly.list() for poly in args[:-1])))]))

def apply_polyvec(f, *args):
    return vector([apply_poly(f, *polys, args[-1]) for polys in zip(*args[:-1])])

## MODULAR REDUCTION

def centered_mod(r, a):
    r_cmod_a = r % a
    if r_cmod_a > a//2:
        r_cmod_a -= a
    return r_cmod_a

## SAMPLING AND HASHING

def sample_matrix_unif(rows, cols, a, b):
    return Matrix([[Rq(sum([random.randint(a, b)*x^i for i in range(n)])) for i_ in range(cols)] for j_ in range(rows)])

def sample_vector_unif(dim, a, b):
    return vector([Rq(sum([random.randint(a, b)*x^i for i in range(n)])) for i_ in range(dim)])

def sample_matrix_hash(rows, cols, a, b, r):
    chash = hashlib.sha256()
    chash.update(r)
    random.seed(int(chash.digest().hex(), 16))
    return Matrix([[Rq(sum([random.randint(a, b)*x^i for i in range(n)])) for i_ in range(cols)] for j_ in range(rows)])

def sample_vector_hash(dim, a, b, pk, tau):
    chash = hashlib.sha256()
    chash.update(bytes(coeff_byte for poly in pk for coeff in poly for coeff_byte in int(coeff).to_bytes(8, byteorder="big")))
    chash.update(int(tau).to_bytes(8, byteorder="big"))
    random.seed(int(chash.digest().hex(), 16))
    return vector([Rq(sum([random.randint(a, b)*x^i for i in range(n)])) for i_ in range(dim)])
    
def sample_challenge_hash(msg, polyvec, pk):
    chash = hashlib.sha256()
    chash.update(bytes(msg, "utf-8"))
    chash.update(bytes(coeff_byte for poly in polyvec for coeff in poly for coeff_byte in int(coeff).to_bytes(8, byteorder="big")))
    chash.update(bytes(coeff_byte for poly in pk for coeff in poly for coeff_byte in int(coeff).to_bytes(8, byteorder="big")))
    random.seed(int(chash.digest().hex(), 16))
    c_list = random.sample([0]*(n-n1c) + random.choices([-1,1], k=n1c), n)
    return Rq(sum([c_list[i]*x^i for i in range(n)]))

## SUPPORTING ALGORITHMS

def power2round(r, u):
    r = r % p
    r0 = centered_mod(r, 2^u)
    return (r - r0)//(2^u)

def decompose(r, alpha):
    r = r % p
    r0 = centered_mod(r, alpha)
    if r - r0 == p - 1:
        r1 = 0
        r0 = r0 - 1
    else:
        r1 = (r - r0)//alpha
        return r1, r0

def highbits(r, alpha):
    r1, r0 = decompose(r, alpha)
    return r1

def lowbits(r, alpha):
    r1, r0 = decompose(r, alpha)
    return r0

def makehint(z, r, alpha):
    r1 = highbits(r, alpha)
    v1 = highbits(r+z, alpha)
    return r1 != v1

def usehint(h, r, alpha):
    m = (p-1)//alpha
    r1, r0 = decompose(r, alpha)
    if h == 1 and r0 > 0:
        return (r1 + 1) % m
    if h == 1 and r0 <= 0:
        return (r1 - 1) % m
    return r1

def power2round_polyvec(r, u):
    return apply_polyvec(power2round, r, u)

def decompose_polyvec(r, alpha):
    return apply_polyvec(decompose, r, alpha)

def highbits_polyvec(r, alpha):
    return apply_polyvec(highbits, r, alpha)

def lowbits_polyvec(r, alpha):
    return apply_polyvec(lowbits, r, alpha)

def makehint_polyvec(z, r, alpha):
    return apply_polyvec(makehint, z, r, alpha)

def usehint_polyvec(h, r, alpha):
    return apply_polyvec(usehint, h, r, alpha)

## NORMS

def infinity_norm_mod(r):
    return abs(centered_mod(lift(r), p))

def infinity_norm_poly(poly):
    return max(infinity_norm_mod(r) for r in poly)

def infinity_norm_polyvec(polyvec):
    return max(infinity_norm_poly(poly) for poly in polyvec)

class blDilithiumQROM_proof:
    @staticmethod
    def identity_keygen(mat):
        K = bytes(random.getrandbits(8) for i in range(seedbytes))
        s1 = sample_vector_unif(l, -eta, eta)
        s2 = sample_vector_unif(k, -eta, eta)
        t = mat*s1 + s2
        t1 = power2round_polyvec(t, d)
        t0 = t - (t1 * (2^d))
        pk = (t1, t0)
        sk = (K, s1, s2)
        return (pk, sk)
    @staticmethod
    def blind_pk(pk, epoch, mat):
        t1, t0 = pk
        s1_p = sample_vector_hash(l, -eta, eta, t1, epoch)
        s2_p = sample_vector_hash(k, -eta, eta, t1, epoch)
        t_p = mat*s1_p + s2_p
        t = t1*(2^d) + t0
        t_bl = t + t_p
        t1_bl = power2round_polyvec(t_bl, d)
        t0_bl = t_bl - (t1_bl * (2^d))
        return (t1_bl, t0_bl)
    @staticmethod
    def sign(sk, pk, M, epoch, mat):
        t1, t0 = pk
        K, s1, s2 = sk
        s1_p = sample_vector_hash(l, -eta, eta, t1, epoch)
        s2_p = sample_vector_hash(k, -eta, eta, t1, epoch)
        s1_bl = s1 + s1_p
        s2_bl = s2 + s2_p
        t_bl = mat*s1_bl + s2_bl
        t1_bl = power2round_polyvec(t_bl, d)
        t0_bl = t_bl - (t1_bl * (2^d))
        kappa = 0
        z, h = "perp", "perp"
        while z == "perp" and h == "perp":
            kappa = kappa + 1
            y = sample_vector_unif(l, -gamma2+1, gamma2-1)
            w = mat*y
            w1 = highbits_polyvec(w, 2*gamma1)
            c = sample_challenge_hash(M, w1, t1_bl)
            z = y + c*s1_bl
            if infinity_norm_polyvec(z) >= gamma2 - beta or infinity_norm_polyvec(lowbits_polyvec(w-c*s2_bl, 2*gamma1)) >= gamma1 - beta:
                z, h = "perp", "perp"
            else:
                h = makehint_polyvec(-1*c*t0_bl, w - c*s2_bl + c*t0_bl, 2*gamma1)
        return (z, h, c)
    @staticmethod
    def verify(pk_bl, M, sigma, mat):
        t1_bl, t0_bl = pk_bl
        z, h, c = sigma
        w1_p = usehint_polyvec(h, mat*z - c*t1_bl*(2^d), 2*gamma1)
        return (infinity_norm_polyvec(z) < gamma2 - beta and c == sample_challenge_hash(M, w1_p, t1_bl))

class blDilithiumQROM:
    @staticmethod
    def identity_keygen(mat):
        K = bytes(random.getrandbits(8) for i in range(seedbytes))
        s1 = sample_vector_unif(l, -eta, eta)
        s2 = sample_vector_unif(k, -eta, eta)
        t = mat*s1 + s2
        t1 = power2round_polyvec(t, d-1)
        t0 = t - (apply_polyvec(lambda a,b: a//b, t1, 2) * (2^d))
        pk = t1
        sk = (K, s1, s2, t0)
        return (pk, sk)
    @staticmethod
    def blind_pk(pk, epoch, mat):
        t1 = pk
        s1_p = sample_vector_hash(l, -eta, eta, t1, epoch)
        s2_p = sample_vector_hash(k, -eta, eta, t1, epoch)
        t_p = mat*s1_p + s2_p
        t1_p = power2round_polyvec(t_p, d-1)
        t1_bl = t1 + t1_p
        # mod step needed?
        return t1_bl
    @staticmethod
    def sign(sk, pk, M, epoch, mat):
        t1 = pk
        K, s1, s2, t0 = sk
        s1_p = sample_vector_hash(l, -eta, eta, t1, epoch)
        s2_p = sample_vector_hash(k, -eta, eta, t1, epoch)
        s1_bl = s1 + s1_p
        s2_bl = s2 + s2_p
        t_bl = mat*s1_bl + s2_bl
        t1_bl = power2round_polyvec(t_bl, d)
        t0_bl = t_bl - (t1_bl * (2^d))
        t1_blp = t1 + power2round_polyvec(mat*s1_p + s2_p, d-1)
        # mod step needed?
        kappa = 0
        z, h = "perp", "perp"
        while z == "perp" and h == "perp":
            kappa = kappa + 1
            y = sample_vector_unif(l, -gamma2+1, gamma2-1)
            w = mat*y
            w1 = highbits_polyvec(w, 2*gamma1)
            c = sample_challenge_hash(M, w1, t1_blp)
            z = y + c*s1_bl
            if infinity_norm_polyvec(z) >= gamma2 - beta or infinity_norm_polyvec(lowbits_polyvec(w-c*s2_bl, 2*gamma1)) >= gamma1 - beta:
                z, h = "perp", "perp"
            else:
                h = makehint_polyvec(-1*c*t0_bl, w - c*s2_bl + c*t0_bl, 2*gamma1)
        return (z, h, c)
    @staticmethod
    def verify(pk_bl, M, sigma, mat):
        t1_bl = pk_bl
        z, h, c = sigma
        w1_p = usehint_polyvec(h, mat*z - c*t1_bl*(2^(d-1)), 2*gamma1)
        return (infinity_norm_polyvec(z) < gamma2 - beta and c == sample_challenge_hash(M, w1_p, t1_bl))

class DilithiumQROM:
    @staticmethod
    def keygen():
        K = bytes(random.getrandbits(8) for i in range(seedbytes))
        rho = bytes(random.getrandbits(8) for i in range(8))
        A = sample_matrix_hash(k, l, 0, p-1, rho)
        s1 = sample_vector_unif(l, -eta, eta)
        s2 = sample_vector_unif(k, -eta, eta)
        t = A*s1 + s2
        t1 = power2round_polyvec(t, d)
        t0 = t - (t1 * (2^d))
        pk = (rho, t1)
        sk = (K, rho, s1, s2, t0)
        return (pk, sk)
    @staticmethod
    def sign(sk, M):
        K, rho, s1, s2, t0 = sk
        A = sample_matrix_hash(k, l, 0, p-1, rho)
        kappa = 0
        z, h = "perp", "perp"
        while z == "perp" and h == "perp":
            kappa = kappa + 1
            y = sample_vector_unif(l, -gamma2+1, gamma2-1)
            w = A*y
            w1 = highbits_polyvec(w, 2*gamma1)
            c = sample_challenge_hash(M, w1, w1)
            z = y + c*s1
            if infinity_norm_polyvec(z) >= gamma2 - (beta//2) or infinity_norm_polyvec(lowbits_polyvec(w-c*s2, 2*gamma1)) >= gamma1 - (beta//2):
                z, h = "perp", "perp"
            else:
                h = makehint_polyvec(-1*c*t0, w - c*s2 + c*t0, 2*gamma1)
        return (z, h, c)
    @staticmethod
    def verify(pk, M, sigma):
        rho, t1 = pk
        A = sample_matrix_hash(k, l, 0, p-1, rho)
        z, h, c = sigma
        w1_p = usehint_polyvec(h, A*z - c*t1*(2^d), 2*gamma1)
        return (infinity_norm_polyvec(z) < gamma2 - (beta//2) and c == sample_challenge_hash(M, w1_p, w1_p))

## CHECK blDilithiumQROM_proof SIGNATURE SCHEME CORRECTNESS

## A = sample_matrix_unif(k, l, 0, p-1)
## pk, sk = blDilithiumQROM_proof.identity_keygen(A)

## check valid signatures verify
## print("vld")
## for i in range(5):
##     message = "".join(random.choices(string.ascii_letters, k=100))
##     ep = random.randint(0, 2^64)
## 
##     pk_bl = blDilithiumQROM_proof.blind_pk(pk, ep, A)
##     sig = blDilithiumQROM_proof.sign(sk, pk, message, ep, A);
##     res = blDilithiumQROM_proof.verify(pk_bl, message, sig, A)
##     if (not res):
##         print(i)

## print("nvld")
## for i in range(5):
##     message = "".join(random.choices(string.ascii_letters, k=100))
##     ep = random.randint(0, 2^64)
## 
##     pk_bl = blDilithiumQROM_proof.blind_pk(pk, ep, A)
##     sig = blDilithiumQROM_proof.sign(sk, pk, message, ep, A);
##     res = blDilithiumQROM_proof.verify(blDilithiumQROM_proof.blind_pk(pk, random.randint(0, 2^64), A), message, sig, A)
##     if (res):
##         print(i)

## CHECK blDilithiumQROM SIGNATURE SCHEME CORRECTNESS
## A = sample_matrix_unif(k, l, 0, p-1)
## pk, sk = blDilithiumQROM.identity_keygen(A)

## check valid signatures verify
## print("vld")
## for i in range(5):
##     message = "".join(random.choices(string.ascii_letters, k=100))
##     ep = random.randint(0, 2^64)
## 
##     pk_bl = blDilithiumQROM.blind_pk(pk, ep, A)
##     sig = blDilithiumQROM.sign(sk, pk, message, ep, A);
##     res = blDilithiumQROM.verify(pk_bl, message, sig, A)
##     if (not res):
##         print(i)

## print("nvld")
## for i in range(5):
##     message = "".join(random.choices(string.ascii_letters, k=100))
##     ep = random.randint(0, 2^64)
## 
##     pk_bl = blDilithiumQROM.blind_pk(pk, ep, A)
##     sig = blDilithiumQROM.sign(sk, pk, message, ep, A);
##     res = blDilithiumQROM.verify(blDilithiumQROM.blind_pk(pk, random.randint(0, 2^64), A), message, sig, A)
##     if (res):
##         print(i)


## CHECK DilithiumQROM SIGNATURE SCHEME CORRECTNESS
## pk, sk = DilithiumQROM.keygen()

## check valid signatures verify
## print("vld")
## for i in range(5):
##     message = "".join(random.choices(string.ascii_letters, k=100))
## 
##     sig = DilithiumQROM.sign(sk, message);
##     res = DilithiumQROM.verify(pk, message, sig)
##     if (not res):
##         print(i)

## print("nvld")
## for i in range(5):
##     message = "".join(random.choices(string.ascii_letters, k=100))
## 
##     sig = DilithiumQROM.sign(sk, message);
##     res = DilithiumQROM.verify(DilithiumQROM.keygen()[0], message, sig)
##     if (res):
##         print(i)

# Number of iterations in the test
rep = 1

print("DilithiumQROM running...")
with open("DilithiumQROM_timings.csv", "w") as file:
    csv_writer = csv.writer(file)
    for i in range(rep):
        time_keygen = 0
        time_blindpk = 0
        time_sign = 0
        time_verify = 0
        message = "".join(random.choices(string.ascii_letters, k=100))

        start = time.time()
        pk, sk = DilithiumQROM.keygen()
        end = time.time()
        time_keygen += (end - start)
        start = time.time()
        sig = DilithiumQROM.sign(sk, message);
        end = time.time()
        time_sign += (end - start)
        start = time.time()
        DilithiumQROM.verify(pk, message, sig)
        end = time.time()
        time_verify += (end - start)
        csv_writer.writerow([time_keygen, time_sign, time_verify])

print("blDilithiumQROM runnning...")
with open("blDilithiumQROM_timings.csv", "w") as file:
    csv_writer = csv.writer(file)
    A = sample_matrix_unif(k, l, 0, p-1)
    for i in range(rep):
        time_keygen = 0
        time_blindpk = 0
        time_sign = 0
        time_verify = 0
        message = "".join(random.choices(string.ascii_letters, k=100))
        ep = random.randint(0, 2^64)

        start = time.time()
        pk, sk = blDilithiumQROM.identity_keygen(A)
        end = time.time()
        time_keygen += (end - start)
        start = time.time()
        pk_bl = blDilithiumQROM.blind_pk(pk, ep, A)
        end = time.time()
        time_blindpk += (end - start)
        start = time.time()
        sig = blDilithiumQROM.sign(sk, pk, message, ep, A);
        end = time.time()
        time_sign += (end - start)
        start = time.time()
        res = blDilithiumQROM.verify(pk_bl, message, sig, A)
        end = time.time()
        time_verify += (end - start)
        csv_writer.writerow([time_keygen, time_blindpk, time_sign, time_verify])

print("DONE!")
