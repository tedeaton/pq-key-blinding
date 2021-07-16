import random
import hashlib
import string

## PARAMETERS

p = 35184372067549
n = 512
k = 4
l = 4
d = 8
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

def apply_poly_coeff_wise(poly, f, *args):
    if len(args) == 0:
        return Rq(sum([f(lift(coeff))*x^i for i, coeff in enumerate(poly.list())]))
    else:
        return Rq(sum([f(lift(coeff), *args)*x^i for i, coeff in enumerate(poly.list())]))

def apply_vector_element_wise(polyvec, f, *args):
    return vector([apply_poly_coeff_wise(poly, f, *args) for poly in polyvec])

def apply_matrix_element_wise(polymat, f, *args):
    return Matrix([apply_vector_element_wise(polyvec, f, *args) for polyvec in polymat])

def apply_poly_binary_wise(poly1, poly2, f, *args):
    if len(args) == 0:
        return Rq(sum([f(lift(coeff1), lift(coeff2))*x^i for i, (coeff1, coeff2) in enumerate(zip(poly1.list(), poly2.list()))]))
    else:
        return Rq(sum([f(lift(coeff1), lift(coeff2), *args)*x^i for i, (coeff1, coeff2) in enumerate(zip(poly1.list(), poly2.list()))]))
    
def apply_vector_binary_wise(polyvec1, polyvec2, f, *args):
    return vector([apply_poly_binary_wise(poly1, poly2, f, *args) for poly1, poly2 in zip(polyvec1, polyvec2)])

def apply_matrix_binary_wise(polymat1, polymat2, f, *args):
    return Matrix([apply_vector_binary_wise(polyvec1, polyvec2, f, *args) for polyvec1, polyvec2 in zip(polymat1, polymat2)])

## MODULAR REDUCTION

def centered_mod(r, a):
    r_cmod_a = r % a
    if r_cmod_a > a//2:
        r_cmod_a -= a
    return r_cmod_a

## SAMPLING AND HASHING

def sample_public_matrix(rows, cols):
    return Matrix([[Rq(P.random_element(degree=n)) for i in range(cols)] for j in range(rows)])

def sample_secret_vector(dim):
    return vector([Rq(sum([random.randint(-eta, eta)*x^i for i in range(n)])) for j in range(dim)])

def sample_blinded_secret_vector(dim, pk, tau):
    chash = hashlib.sha256()
    chash.update(bytes(lift(coeff)%256 for poly in pk[0] for coeff in poly))
    chash.update(bytes(lift(coeff)%256 for poly in pk[1] for coeff in poly))
    chash.update(tau.to_bytes(8, byteorder="big"))
    random.seed(int(chash.digest().hex(), 16))
    return vector([Rq(sum([random.randint(-eta, eta)*x^i for i in range(n)])) for j in range(dim)])

def sample_sign_vector(dim, key, msg, nonce):
    chash = hashlib.sha256()
    chash.update(key)
    chash.update(bytes(msg, "utf-8"))
    chash.update(bytes(nonce))
    random.seed(int(chash.digest().hex(), 16))
    return vector([Rq(sum([random.randint(-gamma1+1, gamma1)*x^i for i in range(n)])) for j in range(dim)])

def sample_challenge(msg, polyvec, pk):
    chash = hashlib.sha256()
    chash.update(bytes(msg, "utf-8"))
    chash.update(bytes(lift(coeff)%256 for poly in polyvec for coeff in poly))
    chash.update(bytes(lift(coeff)%256 for poly in pk[0] for coeff in poly))
    chash.update(bytes(lift(coeff)%256 for poly in pk[1] for coeff in poly))
    random.seed(int(chash.digest().hex(), 16))
    c_list = random.sample([0]*(n-n1c) + random.choices([-1,1], k=n1c), n)
    return Rq(sum([c_list[i]*x^i for i in range(n)]))

## SUPPORTING ALGORITHMS

def power2roundq(r):
    r = r % p
    r0 = centered_mod(r, 2^d)
    return (r - r0)//(2^d)

def decomposeq(r, alpha):
    r = r % p
    r0 = centered_mod(r, alpha)
    if r - r0 == p - 1:
        r1 = 0
        r0 = r0 - 1
    else:
        r1 = (r - r0)//alpha
        return r1, r0

def highbitsq(r, alpha):
    r1, r0 = decomposeq(r, alpha)
    return r1

def lowbitsq(r, alpha):
    r1, r0 = decomposeq(r, alpha)
    return r0

def makehintq(zb, r, alpha):
    r1 = highbitsq(r, alpha)
    v1 = highbitsq(r+zb, alpha)
    return r1 != v1

def usehintq(hb, r, alpha):
    m = (p-1)//alpha
    r1, r0 = decomposeq(r, alpha)
    if hb == 1 and r0 > 0:
        return (r1 + 1) % m
    if hb == 1 and r0 <= 0:
        return (r1 - 1) % m
    return r1

def power2roundq_polyvec(polyvec):
    return apply_vector_element_wise(polyvec, power2roundq)

def decomposeq_polyvec(polyvec, alpha):
    return apply_vector_element_wise(polyvec, decomposeq, alpha)

def highbitsq_polyvec(polyvec, alpha):
    return apply_vector_element_wise(polyvec, highbitsq, alpha)

def lowbitsq_polyvec(polyvec, alpha):
    return apply_vector_element_wise(polyvec, lowbitsq, alpha)

def makehintq_polyvec(polyvec1, polyvec2, alpha):
    return apply_vector_binary_wise(polyvec1, polyvec2, makehintq, alpha)

def usehintq_polyvec(polyvec1, polyvec2, alpha):
    return apply_vector_binary_wise(polyvec1, polyvec2, usehintq, alpha)

## NORMS

def infinity_norm_mod(r):
    return abs(centered_mod(lift(r),p))

def infinity_norm_poly(poly):
    return max(infinity_norm_mod(r) for r in poly)

def infinity_norm_polyvec(polyvec):
    return max(infinity_norm_poly(poly) for poly in polyvec)

## IDENTITY KEYGEN

def identity_keygen(mat):
    K = bytes(random.getrandbits(8) for i in range(seedbytes))
    s1 = sample_secret_vector(l)
    s2 = sample_secret_vector(k)
    t = mat*s1 + s2
    t1 = power2roundq_polyvec(t)
    t0 = t - (t1 * (2^d))
    pk = (t1, t0)
    sk = (K, s1, s2)
    return (pk, sk)

## BLINDING

def blind_pk(pk, epoch, mat):
    t1, t0 = pk
    s1_p = sample_blinded_secret_vector(l, pk, epoch)
    s2_p = sample_blinded_secret_vector(k, pk, epoch)
    t_p = mat*s1_p + s2_p
    t = t1*(2^d) + t0
    t_bl = t + t_p
    t1_bl = power2roundq_polyvec(t_bl)
    t0_bl = t_bl - (t1_bl * (2^d))
    return (t1_bl, t0_bl)

def blind_sk(pk, sk, epoch, mat):
    t1, t0 = pk
    K, s1, s2 = sk
    s1_p = sample_blinded_secret_vector(l, pk, epoch)
    s2_p = sample_blinded_secret_vector(k, pk, epoch)
    s1_bl = s1 + s1_p
    s2_bl = s2 + s2_p
    t_bl = mat*s1_bl + s2_bl
    t1_bl = power2roundq_polyvec(t_bl)
    t0_bl = t_bl - (t1_bl * (2^d))
    pk_bl = (t1_bl, t0_bl)
    sk_bl = (K, s1_bl, s2_bl)
    return (pk_bl, sk_bl)

## SIGN

def sign(sk, pk, M, mat):
    t1, t0 = pk
    K, s1, s2 = sk
    kappa = 0
    z, h = "perp", "perp"
    while z == "perp" and h == "perp":
        kappa = kappa + 1
        y = sample_sign_vector(l, K, M, kappa)
        w = mat*y
        w1 = highbitsq_polyvec(w, 2*gamma2)
        c = sample_challenge(M, w1, pk)
        z = y + c*s1
        if infinity_norm_polyvec(z) >= gamma1 - beta or infinity_norm_polyvec(lowbitsq_polyvec(w-c*s2, 2*gamma2)) >= gamma2 - beta:
            z, h = "perp", "perp"
        else:
            h = makehintq_polyvec(-1*c*t0, w - c*s2 + c*t0, 2*gamma2)
    return (z, h, c)

## VERIFY

def verify(pk, M, sigma, mat):
    t1, t0 = pk
    z, h, c = sigma
    w1_p = usehintq_polyvec(h, mat*z - c*t1*(2^d), 2*gamma2)
    return (infinity_norm_polyvec(z) < gamma1 - beta and c == sample_challenge(M, w1_p, pk))

## CHECK SIGNATURE SCHEME CORRECTNESS

A = sample_public_matrix(k, l)
pk, sk = identity_keygen(A)

## check valid signatures verify
for i in range(5):
    message = "".join(random.choices(string.ascii_letters, k=100))
    epoch = random.randint(0, 2^64)
    pk_bl = blind_pk(pk, epoch, A)
    pk2_bl, sk_bl = blind_sk(pk, sk, epoch, A)
    sig = sign(sk_bl, pk2_bl, message, A);
    res = verify(pk_bl, message, sig, A)
    if (not res):
        print(i)
        
for i in range(5):
    message = "".join(random.choices(string.ascii_letters, k=100))
    epoch = random.randint(0, 2^64)
    pk_bl = blind_pk(pk, epoch, A)
    pk2_bl, sk_bl = blind_sk(pk, sk, epoch, A)
    sig = sign(sk_bl, pk2_bl, message, A);
    res = verify(blind_pk(pk, random.randint(0, 2^64), A), message, sig, A)
    if (res):
        print(i)

print("DONE!")