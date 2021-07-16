#include "blind.h"

void blind_get_challenges(const unsigned char *hash, uint32_t *challenges_index, uint8_t *challenges_sign){
	unsigned char tmp_hash[SEED_BYTES];
	memcpy(tmp_hash,hash,SEED_BYTES);

	// slow hash function
	for(int i=0; i<HASHES; i++){
		HASH(tmp_hash,SEED_BYTES,tmp_hash);
	}

	// generate pseudorandomness
	EXPAND(tmp_hash,SEED_BYTES,(unsigned char *) challenges_index,sizeof(uint32_t)*ROUNDS);

	// set sign bit and zero out higher order bits
	for(int i=0; i<ROUNDS; i++){
		challenges_sign[i] = (challenges_index[i] >> PK_TREE_DEPTH) & 1;
		challenges_index[i] &= (((uint16_t) 1)<<PK_TREE_DEPTH)-1;
	}
}

void blinded_sign(const unsigned char *sk,const unsigned char *m, uint64_t mlen, unsigned char *sig, uint64_t *sig_len, const unsigned char *blinding_nonce)
{
	init_classgroup();

	// hash the message
	unsigned char m_hash[HASH_BYTES];
	HASH(m,mlen,m_hash);

	// pick random seeds
	unsigned char seeds[SEED_BYTES*ROUNDS];
	RAND_bytes(seeds,SEED_BYTES*ROUNDS);

	// compute curves
	mpz_t r[ROUNDS];
	uint curves[ROUNDS] = {{{0}}};


	for(int k=0 ; k<ROUNDS; k++){
		private_key priv;

		// sample mod class number and convert to vector
		mpz_init(r[k]);
		sample_mod_cn_with_seed(seeds + k*SEED_BYTES,r[k]);
		mod_cn_2_vec(r[k],priv.e);

		// compute E_o * vec
		public_key out;
		action(&out, &base, &priv);

		// convert to uint64_t's
		fp_dec(&curves[k], &out.A);
	}

	// hash curves
	unsigned char curve_hash[HASH_BYTES];
	HASH((unsigned char *) curves, sizeof(uint[ROUNDS]), curve_hash);

	// compute master hash
	unsigned char in_buf[2*HASH_BYTES], master_hash[HASH_BYTES];
	memcpy(in_buf,m_hash,HASH_BYTES);
	memcpy(in_buf + HASH_BYTES, curve_hash, HASH_BYTES);
	HASH(in_buf,2*HASH_BYTES, master_hash);

	// copy hash to signature
	memcpy(SIG_HASH(sig),master_hash,HASH_BYTES);

	// get challenges
	uint32_t challenges_index[ROUNDS];
	uint8_t challenges_sign[ROUNDS];
	blind_get_challenges(master_hash,challenges_index,challenges_sign);

	// generate seeds
	unsigned char *sk_seeds = malloc(SEED_BYTES*PKS);
	EXPAND(sk,SEED_BYTES,sk_seeds,SEED_BYTES*PKS);

	unsigned char *update_seeds = malloc(SEED_BYTES*PKS);
	EXPAND(blinding_nonce,SEED_BYTES,update_seeds,SEED_BYTES*PKS);

	// generate secrets mod p
	unsigned char *indices = calloc(1,PKS);
	(void) indices;
	mpz_t s[ROUNDS];
	for(int i=0; i<ROUNDS; i++){

		indices[challenges_index[i]] = 1;
		mpz_init(s[i]);

		mpz_t upd;
		mpz_init(upd);

		sample_mod_cn_with_seed(sk_seeds + challenges_index[i]*SEED_BYTES ,s[i]);
		sample_mod_cn_with_seed(update_seeds + challenges_index[i]*SEED_BYTES, upd);

		mpz_add(s[i], s[i], upd);

		if(challenges_sign[i]){
			mpz_mul_si(s[i],s[i],-1);
		}

		mpz_sub(r[i],s[i],r[i]);
		mpz_fdiv_r(r[i],r[i],cn);

		// silly trick to force export to have 33 bytes
		mpz_add(r[i],r[i],cn);

		mpz_export(SIG_RESPONSES(sig) + 33*i, NULL, 1, 1, 1, 0, r[i]);

		mpz_clear(upd);
		mpz_clear(s[i]);
		mpz_clear(r[i]);
	}

	(*sig_len) = SIG_BYTES;

	clear_classgroup();
	free(indices);
	free(sk_seeds);
	free(update_seeds);
}


void blind_public_key(unsigned char *blinded_pk, const unsigned char *pk, const unsigned char *blinding_nonce)
{

	init_classgroup();

	// generate seeds
	unsigned char *update_seeds = malloc(SEED_BYTES*PKS);
	EXPAND(blinding_nonce,SEED_BYTES,update_seeds,SEED_BYTES*PKS);

	uint* updated_curves = (uint*) PK_CURVES(blinded_pk);
	uint* pkcurves = (uint*) PK_CURVES(pk); 

	for (int i = 0; i < PKS; ++i)
	{
		private_key vec;
		sample_from_classgroup_with_seed(update_seeds + i*SEED_BYTES,vec.e);

		public_key E_i;
		fp_enc(&(E_i.A), &pkcurves[i]);

		// compute E_i * vec
		public_key out;
		action(&out, &E_i, &vec);

		// convert to uint64_t's
		fp_dec(&updated_curves[i], &out.A);
	}

	clear_classgroup();
	free(update_seeds);
}


