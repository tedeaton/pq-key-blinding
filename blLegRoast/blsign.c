#include "blsign.h"

void sign(const unsigned char *sk, const unsigned char *bl, const unsigned char *blpk, const unsigned char *msg, uint64_t msglen, unsigned char *sig, uint64_t *sig_len){
    memset(sig, 0, SIG_BYTES);
    unsigned char buffer[5*HASH_BYTES];

    // hash the message
    HASH(msg, msglen, buffer);

    // Prover part 1
    prover_state state = {0};
    prover_part_1(sk, bl, blpk, SIG_MESSAGE1(sig), &state);

    // Verifier part 1
    HASH(SIG_MESSAGE1(sig), MESSAGE1_BYTES, buffer + HASH_BYTES);
    HASH(buffer, 2*HASH_BYTES, buffer + HASH_BYTES);
    unsigned char challenge1[CHALLENGE1_BYTES];
    generate_challenge1(buffer + HASH_BYTES, challenge1);

    // Prover part 2
    prover_part_2(sk, bl, challenge1, SIG_MESSAGE2(sig), &state);

    // Verifier part 2
    HASH(SIG_MESSAGE2(sig), MESSAGE2_BYTES, buffer + 2*HASH_BYTES);
    HASH(buffer + HASH_BYTES, 2*HASH_BYTES, buffer + 2*HASH_BYTES);
    unsigned char challenge2[CHALLENGE2_BYTES];
    generate_challenge2(buffer + 2*HASH_BYTES, challenge2);

    // Prover part 3
    prover_part_3(sk, bl, challenge2, SIG_MESSAGE3(sig), &state);

    // Verifier part 3
    HASH(SIG_MESSAGE3(sig), MESSAGE3_BYTES, buffer + 3*HASH_BYTES);
    HASH(buffer + 2*HASH_BYTES, 2*HASH_BYTES, buffer + 3*HASH_BYTES);
    unsigned char challenge3[CHALLENGE3_BYTES];
    generate_challenge3(buffer + 3*HASH_BYTES, challenge3);
    
    // Prover part 4
    prover_part_4(sk, bl, blpk, challenge1, challenge2, challenge3, SIG_MESSAGE2(sig), SIG_MESSAGE4(sig), &state);

    // Verifier part 4
    HASH(SIG_MESSAGE4(sig), MESSAGE4_BYTES, buffer + 4*HASH_BYTES);
    HASH(buffer + 3*HASH_BYTES, 2*HASH_BYTES, buffer + 4*HASH_BYTES);
    unsigned char challenge4[CHALLENGE4_BYTES];
    generate_challenge4(buffer + 4*HASH_BYTES, challenge4);

    // Prover part 5
    prover_part_5(sk, bl, blpk, challenge1, challenge2, challenge3, challenge4, SIG_MESSAGE5(sig), &state);

}

void blind(unsigned char *bl, const unsigned char *pk, unsigned char *blpk){
    RAND_bytes(bl, SEED_BYTES);
    blindpk(bl, pk, blpk);
}

int verify(const unsigned char *blpk, const unsigned char *msg, uint64_t msglen, const unsigned char *sig){
    // reconstruct challenges
    unsigned char buffer[5*HASH_BYTES];
    
    HASH(msg, msglen, buffer);

    HASH(SIG_MESSAGE1(sig), MESSAGE1_BYTES, buffer + HASH_BYTES);
    HASH(buffer, 2*HASH_BYTES, buffer + HASH_BYTES);
    unsigned char challenge1[CHALLENGE1_BYTES];
    generate_challenge1(buffer + HASH_BYTES, challenge1);

    HASH(SIG_MESSAGE2(sig), MESSAGE2_BYTES, buffer + 2*HASH_BYTES);
    HASH(buffer + HASH_BYTES, 2*HASH_BYTES, buffer + 2*HASH_BYTES);
    unsigned char challenge2[CHALLENGE2_BYTES];
    generate_challenge2(buffer + 2*HASH_BYTES, challenge2);

    HASH(SIG_MESSAGE3(sig), MESSAGE3_BYTES, buffer + 3*HASH_BYTES);
    HASH(buffer + 2*HASH_BYTES, 2*HASH_BYTES, buffer + 3*HASH_BYTES);
    unsigned char challenge3[CHALLENGE3_BYTES];
    generate_challenge3(buffer + 3*HASH_BYTES, challenge3);

    HASH(SIG_MESSAGE4(sig), MESSAGE4_BYTES, buffer + 4*HASH_BYTES);
    HASH(buffer + 3*HASH_BYTES, 2*HASH_BYTES, buffer + 4*HASH_BYTES);
    unsigned char challenge4[CHALLENGE4_BYTES];
    generate_challenge4(buffer + 4*HASH_BYTES, challenge4);

    return check(blpk, SIG_MESSAGE1(sig), challenge1, SIG_MESSAGE2(sig), challenge2, SIG_MESSAGE3(sig), challenge3, SIG_MESSAGE4(sig), challenge4, SIG_MESSAGE5(sig));
}
