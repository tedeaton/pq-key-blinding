#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "csidh.h"
#include "merkletree.h"
#include "stdint.h"
#include "parameters.h"
#include "classgroup.h"
#include "fp.h"
#include "csifish.h"


void blinded_sign(const unsigned char *sk,const unsigned char *m, uint64_t mlen, unsigned char *sig, uint64_t *sig_len, const unsigned char *blinding_nonce);
void blind_public_key(unsigned char *blinded_pk, const unsigned char *pk, const unsigned char *blinding_nonce);

