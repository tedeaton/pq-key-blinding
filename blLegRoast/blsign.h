#ifndef BLSIGN_H
#define BLSIGN_H 

#include <stdlib.h>

//void *aligned_alloc( size_t alignment, size_t size );

#include "parameters.h"
#include "blzkproof.h"

#define SIG_MESSAGE1(sig) (sig)
#define SIG_MESSAGE2(sig) (SIG_MESSAGE1(sig) + MESSAGE1_BYTES)
#define SIG_MESSAGE3(sig) (SIG_MESSAGE2(sig) + MESSAGE2_BYTES)
#define SIG_MESSAGE4(sig) (SIG_MESSAGE3(sig) + MESSAGE3_BYTES)
#define SIG_MESSAGE5(sig) (SIG_MESSAGE4(sig) + MESSAGE4_BYTES)
#define SIG_BYTES (SIG_MESSAGE5(0) + MESSAGE5_BYTES)

void sign(const unsigned char *sk, const unsigned char *bl, const unsigned char *blpk, const unsigned char *msg, uint64_t msglen, unsigned char *sig, uint64_t *sig_len);
void blind(unsigned char *bl, const unsigned char *pk, unsigned char *blpk);
int verify(const unsigned char *pk, const unsigned char *msg, uint64_t msglen, const unsigned char *sig); 

#endif
