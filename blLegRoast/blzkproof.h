#ifndef BLZKPROOF_H
#define BLZKPROOF_H

#include <openssl/rand.h>
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include "parameters.h"
#include "merkletree.h"

#define SK_BYTES SEED_BYTES

// the order of the shares in memory
#define SHARE_K 0
#define SHARE_T (SHARE_K + 1)
#define SHARES_A_VALUES (SHARE_T + 1)
#define SHARES_B_VALUES (SHARES_A_VALUES + 3)
#define SHARES_C_VALUES (SHARES_B_VALUES + 3)
#define SHARES_OUTPUTS (SHARES_C_VALUES + 3)
#define SHARES_R (SHARES_OUTPUTS + 3)
#define SHARES_PER_PARTY (SHARES_R + RESIDUOSITY_SYMBOLS_PER_ROUND)

typedef __int128 int128_t;
typedef unsigned __int128 uint128_t;

typedef struct
{
    unsigned char seed_trees[ROUNDS][SEED_BYTES*(2*PARTIES-1)];
    uint128_t shares[ROUNDS][PARTIES][SHARES_PER_PARTY];
    uint128_t sums[ROUNDS][SHARES_PER_PARTY];
    uint128_t i_indices[ROUNDS*RESIDUOSITY_SYMBOLS_PER_ROUND];
    uint128_t j_indices[ROUNDS*RESIDUOSITY_SYMBOLS_PER_ROUND];
} prover_state;

#define MESSAGE1_COMMITMENT(mes) (mes)
#define MESSAGE1_DELTA_K(mes) (MESSAGE1_COMMITMENT(mes) + HASH_BYTES)
#define MESSAGE1_DELTA_T(mes) (MESSAGE1_DELTA_K(mes) + ROUNDS*sizeof(uint128_t))
#define MESSAGE1_DELTA_TRIPLE(mes) (MESSAGE1_DELTA_T(mes) + ROUNDS*sizeof(uint128_t))
#define MESSAGE1_BYTES (MESSAGE1_DELTA_TRIPLE(0) + ROUNDS*3*PRIME_BYTES)

#define CHALLENGE1_BYTES (ROUNDS*RESIDUOSITY_SYMBOLS_PER_ROUND*sizeof(uint32_t))

#define MESSAGE2_OUTPUT(mes) (mes)
#define MESSAGE2_BYTES (MESSAGE2_OUTPUT(0) + ROUNDS*RESIDUOSITY_SYMBOLS_PER_ROUND*PRIME_BYTES)

#define CHALLENGE2_LAMBDA(ch) (ch)
#define CHALLENGE2_BYTES (CHALLENGE2_LAMBDA(0) + ROUNDS*RESIDUOSITY_SYMBOLS_PER_ROUND*PRIME_BYTES)

#define MESSAGE3_DELTA_Z(mes) (mes)
#define MESSAGE3_BYTES (MESSAGE3_DELTA_Z(0) + ROUNDS*3*PRIME_BYTES)

#define CHALLENGE3_EPSILON(ch) (ch)
#define CHALLENGE3_BYTES (CHALLENGE3_EPSILON(0) + ROUNDS*3*PRIME_BYTES)

#define MESSAGE4_HASH(mes) (mes)
#define MESSAGE4_ALPHA(mes) (MESSAGE4_HASH(mes) + HASH_BYTES)
#define MESSAGE4_BETA(mes) (MESSAGE4_ALPHA(mes) + ROUNDS*3*PRIME_BYTES)
#define MESSAGE4_BYTES (MESSAGE4_BETA(0) + ROUNDS*3*PRIME_BYTES)

#define CHALLENGE4_BYTES (sizeof(uint32_t[ROUNDS]))

#define MESSAGE5_SEEDS(mes) (mes)
#define MESSAGE5_COMMITMENT(mes) (MESSAGE5_SEEDS(mes) + ROUNDS*PARTY_DEPTH*SEED_BYTES)
#define MESSAGE5_BYTES (MESSAGE5_COMMITMENT(0) + ROUNDS*HASH_BYTES)

#define RESPONSE1_BYTES 1
#define RESPONSE2_BYTES 1


void keygen(unsigned char *pk, unsigned char *sk);
void blindpk(const unsigned char *bl, const unsigned char *pk, unsigned char *blpk);
void prover_part_1(const unsigned char *sk, const unsigned char *bl, const unsigned char *blpk, unsigned char *message1, prover_state *state);
void generate_challenge1(const unsigned char *hash, unsigned char *challenge1);
void prover_part_2(const unsigned char *sk, const unsigned char *bl, const unsigned char *challenge1, unsigned char *message2, prover_state *state);
void generate_challenge2(const unsigned char *hash, unsigned char *challenge2);
void prover_part_3(const unsigned char *sk, const unsigned char *bl, const unsigned char *challenge2, unsigned char *message3, prover_state *state);
void generate_challenge3(const unsigned char *hash, unsigned char *challenge3);
void prover_part_4(const unsigned char *sk, const unsigned char *bl, const unsigned char *blpk, const unsigned char *challenge1, const unsigned char *challenge2, const unsigned char *challenge3, const unsigned char *message2, unsigned char *message4, prover_state *state);
void generate_challenge4(const unsigned char *hash, unsigned char *challenge4);
void prover_part_5(const unsigned char *sk, const unsigned char *bl, const unsigned char *blpk, const unsigned char *challenge1, const unsigned char *challenge2, const unsigned char *challenge3, const unsigned char *challenge4, unsigned char *message5, prover_state *state);

int check(const unsigned char *blpk, const unsigned char *message1, const unsigned char *challenge1, const unsigned char *message2, const unsigned char *challenge2, const unsigned char *message3, const unsigned char *challenge3, const unsigned char *message4, const unsigned char *challenge4, const unsigned char *message5);

unsigned char legendre_symbol(uint128_t *a);


#endif
