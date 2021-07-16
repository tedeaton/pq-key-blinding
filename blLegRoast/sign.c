
#include "sign.h"

static inline
uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}
#define TIC printf("\n"); uint64_t cl = rdtsc();
#define TOC(A) printf("%s cycles = %lu \n",#A ,rdtsc() - cl); cl = rdtsc();

void sign(const unsigned char *sk, const unsigned char *pk, const unsigned char *m, uint64_t mlen, unsigned char *sig, uint64_t *sig_len){
	memset(sig,0,SIG_BYTES);
	unsigned char buffer[4*HASH_BYTES];

	// hash the message
	HASH(m,mlen,buffer);

	// phase 1 
	prover_state state = {0};
	commit(sk, pk, SIG_MESSAGE1(sig), &state);

	// phase 2 
	HASH(SIG_MESSAGE1(sig),MESSAGE1_BYTES,buffer + HASH_BYTES);
	HASH(buffer,2*HASH_BYTES,buffer + HASH_BYTES);
	unsigned char challenge1[CHALLENGE1_BYTES];
	generate_challenge1(buffer + HASH_BYTES, challenge1);

	// phase 3
	respond1(sk, pk, challenge1, SIG_MESSAGE2(sig), &state);

	// phase 4
	HASH(SIG_MESSAGE2(sig), MESSAGE2_BYTES, buffer + 2*HASH_BYTES);
	HASH(buffer + HASH_BYTES,2*HASH_BYTES,buffer + 2*HASH_BYTES);
	unsigned char challenge2[CHALLENGE2_BYTES];
	generate_challenge2(buffer + 2*HASH_BYTES, challenge2); 

	// phase 5
	respond2(sk,pk, challenge1, challenge2, SIG_MESSAGE2(sig), SIG_MESSAGE3(sig), &state);

	// phase 6
	HASH(SIG_MESSAGE3(sig), MESSAGE3_BYTES, buffer + 3*HASH_BYTES);
	HASH(buffer + 2*HASH_BYTES,2*HASH_BYTES,buffer + 3*HASH_BYTES);
	unsigned char challenge3[CHALLENGE3_BYTES];
	generate_challenge3(buffer+ 3*HASH_BYTES, challenge3 ); 

	// phase 7
	respond3(sk,pk, challenge1, challenge2, challenge3, SIG_MESSAGE4(sig), &state);
}

int verify(const unsigned char *pk, const unsigned char *m, uint64_t mlen, const unsigned char *sig){
	// reconstruct challenges
	unsigned char buffer[4*HASH_BYTES];

	HASH(m,mlen,buffer);

	HASH(SIG_MESSAGE1(sig),MESSAGE1_BYTES,buffer + HASH_BYTES);
	HASH(buffer,2*HASH_BYTES,buffer + HASH_BYTES);
	unsigned char challenge1[CHALLENGE1_BYTES];
	generate_challenge1(buffer + HASH_BYTES, challenge1);

	HASH(SIG_MESSAGE2(sig), MESSAGE2_BYTES, buffer + 2*HASH_BYTES);
	HASH(buffer + HASH_BYTES,2*HASH_BYTES,buffer + 2*HASH_BYTES);
	unsigned char challenge2[CHALLENGE2_BYTES];
	generate_challenge2(buffer + 2*HASH_BYTES, challenge2); 

	HASH(SIG_MESSAGE3(sig), MESSAGE3_BYTES, buffer + 3*HASH_BYTES);
	HASH(buffer + 2*HASH_BYTES,2*HASH_BYTES,buffer + 3*HASH_BYTES);
	unsigned char challenge3[CHALLENGE3_BYTES];
	generate_challenge3(buffer+ 3*HASH_BYTES, challenge3 ); 

	return check(pk, SIG_MESSAGE1(sig), challenge1, SIG_MESSAGE2(sig), challenge2, SIG_MESSAGE3(sig), challenge3, SIG_MESSAGE4(sig));
}
