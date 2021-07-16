from math import comb, log, sqrt
from decimal import *

# are we trying to find a blinded or unblinded parameter set
blinded = True

# returns the probability that with <trials>
# independent trials with success probability
# <prob>, there are exactly <suc> successes.
def binomial_pdf(trials, prob, suc):
    if suc > trials: return 0
    elif suc < 0: return 0
    with localcontext() as ctx:
        ctx.prec = 1024
        result = Decimal(comb(trials, suc))
        result *= Decimal(prob**suc)
        result *= Decimal((1 - prob)**(trials - suc))
    return result

# returns the probability that with <trials>
# independent trials with success probability
# <prob>, there are at least <suc> successes.
def binomial_cdf(trials, prob, suc):
    with localcontext() as ctx:
        ctx.prec = 1024
        result = Decimal(0)
        for i in range(suc, trials + 1):
            result += binomial_pdf(trials, prob, i)
    return result 


def eq_1_quantity(M, beta, B, N, M_primed):
    with localcontext() as ctx:
        ctx.prec = 1024
        result = Decimal(N**(M - M_primed))
    result += Decimal(binomial_cdf(M, (1 - beta)**B, M_primed))**(-1)
    return result


def is_eq1_satisfied(M, beta, B, N):
    for M_primed in range(0,M+1):
        if eq_1_quantity(M, beta, B, N, M_primed) < Decimal(2**128):
            return False
    return True



def signature_size(N, B, M):
    PRIME_BYTES = 16
    HASH_BYTES = 32
    SEED_BYTES = 16
    

    result = 0
    # message1 bytes
    result += HASH_BYTES # the message 1 commitment
    result += M*PRIME_BYTES # the delta K values
    result += M*PRIME_BYTES # the delta T values
    result += M*3*PRIME_BYTES # the delta triple values

    # message2 bytes
    result += M*B*PRIME_BYTES
    
    # message3 bytes
    result += M*3*PRIME_BYTES # the delta z values

    # message4 bytes
    result += HASH_BYTES # the message 4 commitement
    result += M*3*PRIME_BYTES # the alpha values
    result += M*3*PRIME_BYTES # the beta values

    # message5 bytes
    result += M*log(N, 2)*SEED_BYTES # revealing seeds
    result += M*HASH_BYTES # commitments

    return result



def seek_params(beta):
    for log_N in range(6, 20):
        N = 2**log_N

        min_size = 999999999
        min_params = [0, 0, 0]
        for M in range(100, 10, -1):
            for B in range(10, 100):
                if signature_size(N,B,M) < min_size:
                    if is_eq1_satisfied(M, beta, B, N):
                        min_size = signature_size(N,B,M)
                        min_params = [N, M, B]
                        break
        print("parameter set found:")
        print("N = " + str(min_params[0]) + ", M = " + str(min_params[1]) + ", B = " + str(min_params[2]))
        print("signature size: " + str( min_size))


def is_beta_bound_satisfied(i,L):
    with localcontext() as ctx:
        k = Decimal(254)
        prime = Decimal(2**127 - 1)
        prob = Decimal( 1/k + 1/Decimal(sqrt(prime)) + 2/prime)
        return binomial_cdf(L, prob, i) < Decimal(2**(-390))

def find_beta(L):
    i = L
    while(is_beta_bound_satisfied(i, L) and i > 0):
        i -= 1
    
    return 1.0 - (1.0*i / L)



#seek_params(0.15625)
#print(find_beta(512))
seek_params(0.822)