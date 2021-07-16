#include "zkproof.h"

static inline
uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}
#define TIC printf(""); uint64_t cl = rdtsc();
#define TOC(A) printf("%s cycles = %lu \n",#A ,rdtsc() - cl); cl = rdtsc();

uint128_t m127 = ((((uint128_t) 1) << 127)-1);

inline 
void reduce_mod_p(uint128_t * a){
	while( *a >= m127 ){
		*a -= m127;
	}
}

inline
void add_mod_p(uint128_t *a,uint128_t *b){
	*a += *b;
	if(*a < *b){
		reduce_mod_p(a);
		*a += 2;
	}
}

inline
uint128_t compute_index(uint32_t a){
	uint128_t out = 0;
	EXPAND((unsigned char *)&a, sizeof(a), (unsigned char *) &out, sizeof(uint128_t));
	return out;
}

inline
void compute_indices(const uint32_t *a, uint128_t *indices){
	for (int i = 0; i < ROUNDS*RESIDUOSITY_SYMBOLS_PER_ROUND; ++i)
	{
		indices[i] = compute_index(a[i]);
	}
}

void square_mod_p(uint128_t *out, uint128_t* a){
	reduce_mod_p(a); // TODO: can we remove/avoid this?

	uint128_t lowa  = *((uint64_t *) a);
	uint128_t higha = *((uint64_t *) a + 1);

	*out   = lowa * lowa;
	uint128_t out64  = (lowa * higha) << 1;
	uint128_t out127 = (higha * higha + (out64 >> 64)) << 1 ;

	out64 <<= 64;

	add_mod_p(out,&out127);
	add_mod_p(out,&out64);
}

void mul_add_mod_p(uint128_t* out, uint128_t* a, uint128_t* b){
	reduce_mod_p(a); // TODO: can we remove/avoid this?
	reduce_mod_p(b);

	uint128_t lowa  = (*a) % (((uint128_t) 1) << 64);
	uint128_t lowb  = (*b) % (((uint128_t) 1) << 64);
	uint128_t higha = (*a) >> 64;
	uint128_t highb = (*b) >> 64;

	uint128_t out0   = lowa * lowb;
	uint128_t out64  = (lowa * highb) + (lowb * higha);
	uint128_t out127 = (higha * highb + (out64 >> 64)) << 1 ;

	out64 <<= 64;

	add_mod_p(out,&out0);
	add_mod_p(out,&out127);
	add_mod_p(out,&out64);
	//TODO optimize if necessary ?
}

inline 
void sample_mod_p(const unsigned char *seed, uint128_t *out){
	EXPAND(seed,SEED_BYTES,(unsigned char *) out,sizeof(uint128_t));
	reduce_mod_p(out);
}

void print_uint128(uint128_t n)
{
    if (n == 0) {
      return;
    }

    print_uint128(n/10);
    putchar(n%10+0x30);
}

void printbinary(uint128_t a){
	for (int i = 0; i < 128; ++i){
		printf("%d", a%2);
		a >>= 1;
	}
}

unsigned char legendre_symbol_ct(uint128_t *a){
	uint128_t out = *a;
	uint128_t temp, temp2;

	square_mod_p(&temp,&out);
	out = 0;
	mul_add_mod_p(&out,&temp,a);
	square_mod_p(&temp,&out);
	out = 0;
	mul_add_mod_p(&out,&temp,a);
	square_mod_p(&temp,&out);
	out = 0;
	mul_add_mod_p(&out,&temp,a);
	square_mod_p(&temp,&out);
	out = 0;
	mul_add_mod_p(&out,&temp,a);
	square_mod_p(&temp,&out);
	out = 0;
	mul_add_mod_p(&out,&temp,a);
	temp2 = out;

	for (int i = 0; i < 20; ++i)
	{
		square_mod_p(&temp,&out);
		square_mod_p(&temp,&temp);
		square_mod_p(&temp,&temp);
		square_mod_p(&temp,&temp);
		square_mod_p(&temp,&temp);
		square_mod_p(&temp,&temp);

		out = 0;
		mul_add_mod_p(&out, &temp,&temp2);  
	}

	reduce_mod_p(&out);
	return (-out+1)/2;
}

void test(){
	uint128_t x = ((uint128_t) 1) << 100;
	TIC
	for (int j = 0; j < 20; ++j)
	{
		int sum = 0;
		for (uint128_t i = j*100000; i < 100000+ j*100000; ++i)
		{
			uint128_t a = i+x;
			sum += legendre_symbol_ct(&a);
		}
		printf("sum: %d \n", sum);
		TOC(test);
	}
}

unsigned char power_residue_symbol(uint128_t *a){
	uint128_t out = *a;
	uint128_t temp;
	for (int i = 0; i < 17; ++i)
	{
		// square 7 times and multiply by a
		square_mod_p(&temp,&out);
		square_mod_p(&out,&temp);
		square_mod_p(&temp,&out);
		square_mod_p(&out,&temp);
		square_mod_p(&temp,&out);
		square_mod_p(&out,&temp);
		square_mod_p(&temp,&out);
		out = 0;
		mul_add_mod_p(&out,&temp,a);
	}
	reduce_mod_p(&out);

	static const uint64_t list[2*254] = {1U, 0U, 18446726481523507199U, 9223372036854775807U, 0U, 16777216U, 18446744073709551583U, 9223372036854775807U, 562949953421312U, 0U, 18446744073709551615U, 9223372036317904895U, 1024U, 0U, 18428729675200069631U, 9223372036854775807U, 0U, 17179869184U, 18446744073709518847U, 9223372036854775807U, 576460752303423488U, 0U, 18446744073709551615U, 9223371487098961919U, 1048576U, 0U, 18446744073709551615U, 9223372036854775806U, 0U, 17592186044416U, 18446744073675997183U, 9223372036854775807U, 0U, 32U, 18446744073709551615U, 9222809086901354495U, 1073741824U, 0U, 18446744073709551615U, 9223372036854774783U, 0U, 18014398509481984U, 18446744039349813247U, 9223372036854775807U, 0U, 32768U, 18446744073709551615U, 8646911284551352319U, 1099511627776U, 0U, 18446744073709551615U, 9223372036853727231U, 2U, 0U, 18446708889337462783U, 9223372036854775807U, 0U, 33554432U, 18446744073709551551U, 9223372036854775807U, 1125899906842624U, 0U, 18446744073709551615U, 9223372035781033983U, 2048U, 0U, 18410715276690587647U, 9223372036854775807U, 0U, 34359738368U, 18446744073709486079U, 9223372036854775807U, 1152921504606846976U, 0U, 18446744073709551615U, 9223370937343148031U, 2097152U, 0U, 18446744073709551615U, 9223372036854775805U, 0U, 35184372088832U, 18446744073642442751U, 9223372036854775807U, 0U, 64U, 18446744073709551615U, 9222246136947933183U, 2147483648U, 0U, 18446744073709551615U, 9223372036854773759U, 0U, 36028797018963968U, 18446744004990074879U, 9223372036854775807U, 0U, 65536U, 18446744073709551615U, 8070450532247928831U, 2199023255552U, 0U, 18446744073709551615U, 9223372036852678655U, 4U, 0U, 18446673704965373951U, 9223372036854775807U, 0U, 67108864U, 18446744073709551487U, 9223372036854775807U, 2251799813685248U, 0U, 18446744073709551615U, 9223372034707292159U, 4096U, 0U, 18374686479671623679U, 9223372036854775807U, 0U, 68719476736U, 18446744073709420543U, 9223372036854775807U, 2305843009213693952U, 0U, 18446744073709551615U, 9223369837831520255U, 4194304U, 0U, 18446744073709551615U, 9223372036854775803U, 0U, 70368744177664U, 18446744073575333887U, 9223372036854775807U, 0U, 128U, 18446744073709551615U, 9221120237041090559U, 4294967296U, 0U, 18446744073709551615U, 9223372036854771711U, 0U, 72057594037927936U, 18446743936270598143U, 9223372036854775807U, 0U, 131072U, 18446744073709551615U, 6917529027641081855U, 4398046511104U, 0U, 18446744073709551615U, 9223372036850581503U, 8U, 0U, 18446603336221196287U, 9223372036854775807U, 0U, 134217728U, 18446744073709551359U, 9223372036854775807U, 4503599627370496U, 0U, 18446744073709551615U, 9223372032559808511U, 8192U, 0U, 18302628885633695743U, 9223372036854775807U, 0U, 137438953472U, 18446744073709289471U, 9223372036854775807U, 4611686018427387904U, 0U, 18446744073709551615U, 9223367638808264703U, 8388608U, 0U, 18446744073709551615U, 9223372036854775799U, 0U, 140737488355328U, 18446744073441116159U, 9223372036854775807U, 0U, 256U, 18446744073709551615U, 9218868437227405311U, 8589934592U, 0U, 18446744073709551615U, 9223372036854767615U, 0U, 144115188075855872U, 18446743798831644671U, 9223372036854775807U, 0U, 262144U, 18446744073709551615U, 4611686018427387903U, 8796093022208U, 0U, 18446744073709551615U, 9223372036846387199U, 16U, 0U, 18446462598732840959U, 9223372036854775807U, 0U, 268435456U, 18446744073709551103U, 9223372036854775807U, 9007199254740992U, 0U, 18446744073709551615U, 9223372028264841215U, 16384U, 0U, 18158513697557839871U, 9223372036854775807U, 0U, 274877906944U, 18446744073709027327U, 9223372036854775807U, 9223372036854775808U, 0U, 18446744073709551615U, 9223363240761753599U, 16777216U, 0U, 18446744073709551615U, 9223372036854775791U, 0U, 281474976710656U, 18446744073172680703U, 9223372036854775807U, 0U, 512U, 18446744073709551615U, 9214364837600034815U, 17179869184U, 0U, 18446744073709551615U, 9223372036854759423U, 0U, 288230376151711744U, 18446743523953737727U, 9223372036854775807U, 0U, 524288U, 18446744073709551614U, 9223372036854775807U, 17592186044416U, 0U, 18446744073709551615U, 9223372036837998591U, 32U, 0U, 18446181123756130303U, 9223372036854775807U, 0U, 536870912U, 18446744073709550591U, 9223372036854775807U, 18014398509481984U, 0U, 18446744073709551615U, 9223372019674906623U, 32768U, 0U, 17870283321406128127U, 9223372036854775807U, 0U, 549755813888U, 18446744073708503039U, 9223372036854775807U, 0U, 1U, 18446744073709551615U, 9223354444668731391U, 33554432U, 0U, 18446744073709551615U, 9223372036854775775U, 0U, 562949953421312U, 18446744072635809791U, 9223372036854775807U, 0U, 1024U, 18446744073709551615U, 9205357638345293823U, 34359738368U, 0U, 18446744073709551615U, 9223372036854743039U, 0U, 576460752303423488U, 18446742974197923839U, 9223372036854775807U, 0U, 1048576U, 18446744073709551613U, 9223372036854775807U, 35184372088832U, 0U, 18446744073709551615U, 9223372036821221375U, 64U, 0U, 18445618173802708991U, 9223372036854775807U, 0U, 1073741824U, 18446744073709549567U, 9223372036854775807U, 36028797018963968U, 0U, 18446744073709551615U, 9223372002495037439U, 65536U, 0U, 17293822569102704639U, 9223372036854775807U, 0U, 1099511627776U, 18446744073707454463U, 9223372036854775807U, 0U, 2U, 18446744073709551615U, 9223336852482686975U, 67108864U, 0U, 18446744073709551615U, 9223372036854775743U, 0U, 1125899906842624U, 18446744071562067967U, 9223372036854775807U, 0U, 2048U, 18446744073709551615U, 9187343239835811839U, 68719476736U, 0U, 18446744073709551615U, 9223372036854710271U, 0U, 1152921504606846976U, 18446741874686296063U, 9223372036854775807U, 0U, 2097152U, 18446744073709551611U, 9223372036854775807U, 70368744177664U, 0U, 18446744073709551615U, 9223372036787666943U, 128U, 0U, 18444492273895866367U, 9223372036854775807U, 0U, 2147483648U, 18446744073709547519U, 9223372036854775807U, 72057594037927936U, 0U, 18446744073709551615U, 9223371968135299071U, 131072U, 0U, 16140901064495857663U, 9223372036854775807U, 0U, 2199023255552U, 18446744073705357311U, 9223372036854775807U, 0U, 4U, 18446744073709551615U, 9223301668110598143U, 134217728U, 0U, 18446744073709551615U, 9223372036854775679U, 0U, 2251799813685248U, 18446744069414584319U, 9223372036854775807U, 0U, 4096U, 18446744073709551615U, 9151314442816847871U, 137438953472U, 0U, 18446744073709551615U, 9223372036854644735U, 0U, 2305843009213693952U, 18446739675663040511U, 9223372036854775807U, 0U, 4194304U, 18446744073709551607U, 9223372036854775807U, 140737488355328U, 0U, 18446744073709551615U, 9223372036720558079U, 256U, 0U, 18442240474082181119U, 9223372036854775807U, 0U, 4294967296U, 18446744073709543423U, 9223372036854775807U, 144115188075855872U, 0U, 18446744073709551615U, 9223371899415822335U, 262144U, 0U, 13835058055282163711U, 9223372036854775807U, 0U, 4398046511104U, 18446744073701163007U, 9223372036854775807U, 0U, 8U, 18446744073709551615U, 9223231299366420479U, 268435456U, 0U, 18446744073709551615U, 9223372036854775551U, 0U, 4503599627370496U, 18446744065119617023U, 9223372036854775807U, 0U, 8192U, 18446744073709551615U, 9079256848778919935U, 274877906944U, 0U, 18446744073709551615U, 9223372036854513663U, 0U, 4611686018427387904U, 18446735277616529407U, 9223372036854775807U, 0U, 8388608U, 18446744073709551599U, 9223372036854775807U, 281474976710656U, 0U, 18446744073709551615U, 9223372036586340351U, 512U, 0U, 18437736874454810623U, 9223372036854775807U, 0U, 8589934592U, 18446744073709535231U, 9223372036854775807U, 288230376151711744U, 0U, 18446744073709551615U, 9223371761976868863U, 524288U, 0U, 9223372036854775807U, 9223372036854775807U, 0U, 8796093022208U, 18446744073692774399U, 9223372036854775807U, 0U, 16U, 18446744073709551615U, 9223090561878065151U, 536870912U, 0U, 18446744073709551615U, 9223372036854775295U, 0U, 9007199254740992U, 18446744056529682431U, 9223372036854775807U, 0U, 16384U, 18446744073709551615U, 8935141660703064063U, 549755813888U, 0U, 18446744073709551615U, 9223372036854251519U};
	
	reduce_mod_p(&out);
	for (int i = 0; i < 254; ++i)
	{
		if( memcmp((unsigned char *) &out, (unsigned char *) (list + 2*i), sizeof(uint128_t)) == 0 ){
			return i;
		}
	}
	// oops
	return 0;
}

void keygen(unsigned char *pk, unsigned char *sk){
	RAND_bytes(sk,SEED_BYTES);
	
	uint128_t key;
	sample_mod_p(sk, &key);

	memset(pk,0,PK_BYTES);
	#ifdef LEGENDRE
		for (uint32_t i = 0; i < PK_BYTES*8 ; ++i)
		{
			uint128_t temp = compute_index(i);
			add_mod_p(&temp,&key);
			pk[i/8] |= legendre_symbol_ct(&temp) << (i%8);
		}
	#else
		for (uint32_t i = 0; i < PK_BYTES ; ++i)
		{
			uint128_t temp = compute_index(i);
			add_mod_p(&temp,&key);
			pk[i] = power_residue_symbol(&temp);
		}
	#endif
}

void commit(const unsigned char *sk, const unsigned char *pk, unsigned char *message1, prover_state *state){
	memset(message1,0,MESSAGE1_BYTES);

	// generate key
	uint128_t key;
	sample_mod_p(sk, &key);

	unsigned char commitments[ROUNDS*PARTIES*HASH_BYTES + ROUNDS*RESIDUOSITY_SYMBOLS_PER_ROUND];

	for (int round = 0; round < ROUNDS; round++)
	{
		// pick root seed
		RAND_bytes(state->seed_trees[round],SEED_BYTES);
		// generate seeds
		generate_seed_tree(state->seed_trees[round]);

		// generate the commitments and the shares
		for (int i = 0; i < PARTIES; ++i)
		{
			// commit to seed
			HASH(&state->seed_trees[round][(PARTIES-1+i)*SEED_BYTES],SEED_BYTES,commitments + round*PARTIES*HASH_BYTES + i*HASH_BYTES);

			// generate shares from seed
			EXPAND(&state->seed_trees[round][(PARTIES-1+i)*SEED_BYTES],SEED_BYTES,(unsigned char *) state->shares[round][i],sizeof(state->shares[round][i]));

			// add the shares to the sums
			for (int j = 0; j < SHARES_PER_PARTY; ++j)
			{
				add_mod_p(&state->sums[round][j], &state->shares[round][i][j]);
			}
		}

		// reduce sums mod p
		for (int i = 0; i < SHARES_PER_PARTY; ++i)
		{
			reduce_mod_p(&state->sums[round][i]);
		}

		// compute legendre symbols of R_i
		for (int i = 0; i < RESIDUOSITY_SYMBOLS_PER_ROUND; ++i)
		{
			#ifdef LEGENDRE
				commitments[ROUNDS*PARTIES*HASH_BYTES + round*RESIDUOSITY_SYMBOLS_PER_ROUND + i] = legendre_symbol_ct(&state->sums[round][SHARES_R + i]); 
			#else
				commitments[ROUNDS*PARTIES*HASH_BYTES + round*RESIDUOSITY_SYMBOLS_PER_ROUND + i] = power_residue_symbol(&state->sums[round][SHARES_R + i]);
			#endif
		}

		// compute Delta K and add to share 0
		uint128_t *delta_k = (uint128_t*) MESSAGE1_DELTA_K(message1) + round;
		*delta_k = m127 - state->sums[round][SHARE_K];
		add_mod_p(delta_k , &key);
		reduce_mod_p(delta_k);
		add_mod_p(&(state->shares[round][0][SHARE_K]),delta_k);
		
		// compute Delta Triple and add to share 0
		uint128_t *delta_triple = (uint128_t*) MESSAGE1_DELTA_TRIPLE(message1) + round ;
		*delta_triple = m127 - state->sums[round][SHARES_TRIPLE+2];
		mul_add_mod_p(delta_triple,&(state->sums[round][SHARES_TRIPLE]),&(state->sums[round][SHARES_TRIPLE+1]));
		reduce_mod_p(delta_triple);
		add_mod_p(&(state->shares[round][0][SHARES_TRIPLE+2]),delta_triple);
	}

	HASH(commitments,sizeof(commitments),MESSAGE1_COMMITMENT(message1));
}

void generate_challenge1(const unsigned char *hash, unsigned char *challenge1 ){
	uint32_t *indices = (uint32_t *) challenge1;
	EXPAND(hash, HASH_BYTES, (unsigned char *) indices,sizeof(uint32_t[ROUNDS*RESIDUOSITY_SYMBOLS_PER_ROUND]));
	for (int i = 0; i < ROUNDS*RESIDUOSITY_SYMBOLS_PER_ROUND; ++i)
	{
		#ifdef LEGENDRE
			indices[i] &= (( (uint32_t) 1 )<<PK_DEPTH)-1;
		#else
			indices[i] &= (( (uint32_t) 1 )<<(PK_DEPTH-3))-1;
		#endif
	}
}

void respond1(const unsigned char *sk, const unsigned char *pk, const unsigned char * challenge1, unsigned char *message2, prover_state *state){
	compute_indices((uint32_t *) challenge1, state->indices);

	memset(message2,0,MESSAGE2_BYTES);

	// generate key
	uint128_t key;
	sample_mod_p(sk, &key);

	uint128_t *output = (uint128_t *) MESSAGE2_OUTPUT(message2);

	// Compute output
	for (int round = 0; round < ROUNDS; ++round)
	{
		for (int i = 0; i < RESIDUOSITY_SYMBOLS_PER_ROUND; ++i)
		{	
			uint128_t key_plus_I = state->indices[round*RESIDUOSITY_SYMBOLS_PER_ROUND + i];
			add_mod_p(&key_plus_I, &key);
			mul_add_mod_p(&(output[round*RESIDUOSITY_SYMBOLS_PER_ROUND + i]),&key_plus_I,&(state->sums[round][SHARES_R + i]));
			reduce_mod_p(&(output[round*RESIDUOSITY_SYMBOLS_PER_ROUND + i]));
		}
	} 
}

void generate_challenge2(const unsigned char *hash, unsigned char *challenge2){
	EXPAND(hash,HASH_BYTES,challenge2,CHALLENGE2_BYTES);
	//memset(challenge2,0,CHALLENGE2_BYTES);
}

void respond2(const unsigned char *sk, const unsigned char *pk, const unsigned char *challenge1 , const unsigned char *challenge2, unsigned char *message2, unsigned char *message3, prover_state *state){
	memset(message3,0,MESSAGE3_BYTES);

	// generate key
	uint128_t key;
	sample_mod_p(sk, &key );

	uint128_t openings[ROUNDS][PARTIES][3] = {0};

	const uint32_t *indices = (uint32_t *) challenge1;

	for (int round = 0; round < ROUNDS; ++round)
	{
		uint128_t *epsilon, *alpha, *beta;

		epsilon = (uint128_t*) CHALLENGE2_EPSILON(challenge2) + round;
		alpha = (uint128_t*) MESSAGE3_ALPHA(message3) + round;
		beta = (uint128_t*) MESSAGE3_BETA(message3) + round;

		// compute alpha in the clear
		*alpha = 0;
		mul_add_mod_p(alpha,epsilon,&key);
		add_mod_p(alpha,&state->sums[round][SHARES_TRIPLE]);
		reduce_mod_p(alpha);

		// import lambda and compute beta in the clear
		uint128_t *lambda = (uint128_t*) CHALLENGE2_LAMBDA(challenge2) + round*RESIDUOSITY_SYMBOLS_PER_ROUND;

		*beta = state->sums[round][SHARES_TRIPLE+1];
		for (int i = 0; i < RESIDUOSITY_SYMBOLS_PER_ROUND; ++i)
		{
			mul_add_mod_p(beta,lambda + i, &state->sums[round][SHARES_R + i]);
		}
		reduce_mod_p(beta);

		// computes shares of alpha, beta and v
		for (int i = 0; i < PARTIES; ++i)
		{
			// compute share of alpha
			mul_add_mod_p(&openings[round][i][0],epsilon,&state->shares[round][i][SHARE_K]);
			add_mod_p(&openings[round][i][0],&state->shares[round][i][SHARES_TRIPLE]);
			reduce_mod_p(&openings[round][i][0]);

			// compute share of beta and z
			uint128_t z_share = 0;
			openings[round][i][1] = state->shares[round][i][SHARES_TRIPLE+1];
			for (int j = 0; j < RESIDUOSITY_SYMBOLS_PER_ROUND; ++j)
			{
				// share of beta
				mul_add_mod_p(&openings[round][i][1],&state->shares[round][i][SHARES_R+j],lambda + j);
				reduce_mod_p(&openings[round][i][1]);

				// share of z
				uint128_t temp2 = 0;
				mul_add_mod_p(&temp2,&state->shares[round][i][SHARES_R + j],lambda + j); // this multiplication is done earlier, reuse result
				reduce_mod_p(&temp2);
				temp2 = m127 - temp2;
				uint128_t index = state->indices[round*RESIDUOSITY_SYMBOLS_PER_ROUND + j];
				mul_add_mod_p(&z_share,&temp2,&index); 

				if(i==0){
					mul_add_mod_p(&z_share, lambda+j , ((uint128_t *)MESSAGE2_OUTPUT(message2)) + round*RESIDUOSITY_SYMBOLS_PER_ROUND + j  );
				}
			}

			// compute sharing of v
			openings[round][i][2] = state->shares[round][i][SHARES_TRIPLE+ 2];
			if(i==0){
				mul_add_mod_p(&openings[round][i][2],alpha,beta);
			}
			reduce_mod_p(&openings[round][i][2]);
			openings[round][i][2] = m127 - openings[round][i][2];
			mul_add_mod_p(&openings[round][i][2],alpha,&state->shares[round][i][SHARES_TRIPLE+1]);
			mul_add_mod_p(&openings[round][i][2],beta ,&state->shares[round][i][SHARES_TRIPLE]);
			mul_add_mod_p(&openings[round][i][2],epsilon,&z_share);
			reduce_mod_p(&openings[round][i][2]);
		}
	}

	HASH((unsigned char *)openings,sizeof(openings),MESSAGE3_HASH(message3));
}

void generate_challenge3(const unsigned char *hash, unsigned char *challenge3){
	uint32_t *unopened_party = (uint32_t *) challenge3;
	EXPAND(hash, HASH_BYTES, challenge3,sizeof(uint32_t[ROUNDS]));
	for (int i = 0; i < ROUNDS; ++i)
	{
		unopened_party[i] &= (( (uint32_t) 1 )<<PARTY_DEPTH)-1;
		//unopened_party[i] = 0;
	}
}

void respond3(const unsigned char *sk, const unsigned char *pk, const unsigned char *challenge1, const unsigned char *challenge2, const unsigned char *challenge3, unsigned char *message4, prover_state *state){
	uint32_t *unopened_party = (uint32_t *) challenge3;
	for (int round = 0; round < ROUNDS; ++round)
	{
		release_seeds(state->seed_trees[round],unopened_party[round],MESSAGE4_SEEDS(message4) + round*PARTY_DEPTH*SEED_BYTES);

		HASH(&state->seed_trees[round][(PARTIES-1+unopened_party[round])*SEED_BYTES], SEED_BYTES, MESSAGE4_COMMITMENT(message4) + round*HASH_BYTES);
	}
}

int getBit(const unsigned char *pk, uint32_t bit){
	#ifdef LEGENDRE
		return( (pk[bit/8] & (1 << (bit%8)) ) > 0); 
	#else
		return pk[bit];
	#endif
}

int check(const unsigned char *pk, const unsigned char *message1, const unsigned char *challenge1, const unsigned char *message2, const unsigned char *challenge2, const unsigned char *message3, const unsigned char *challenge3, const unsigned char *message4){
	prover_state state={0};
	compute_indices((uint32_t *) challenge1, state.indices);
	
	uint32_t *unopened_party = (uint32_t *) challenge3;

	// check first commitment: seeds + jacobi syymbols
	unsigned char buffer1[ROUNDS*PARTIES*HASH_BYTES + ROUNDS*RESIDUOSITY_SYMBOLS_PER_ROUND];
	for (int round = 0; round < ROUNDS; ++round)
	{
		fill_down(state.seed_trees[round],unopened_party[round], MESSAGE4_SEEDS(message4) + round*PARTY_DEPTH*SEED_BYTES);

		//copy the commitment of the unopened value
		memcpy(buffer1 + round*PARTIES*HASH_BYTES + unopened_party[round] * HASH_BYTES , MESSAGE4_COMMITMENT(message4) + round*HASH_BYTES , HASH_BYTES);

		// generate the commitments to the remaining seeds
		for (int i = 0; i < PARTIES; ++i)
		{
			if(i == unopened_party[round]){
				continue;
			}

			// commit to seed
			HASH(&(state.seed_trees[round][(PARTIES-1+i)*SEED_BYTES]),SEED_BYTES,buffer1 + round*PARTIES*HASH_BYTES + i*HASH_BYTES);

			// generate shares from seed (and Delta's)
			EXPAND(&state.seed_trees[round][(PARTIES-1+i)*SEED_BYTES],SEED_BYTES,(unsigned char *) state.shares[round][i],sizeof(state.shares[round][i]));
			if(i == 0){
				add_mod_p(&state.shares[round][i][SHARE_K], ((uint128_t *) MESSAGE1_DELTA_K(message1) ) + round  );
				add_mod_p(&state.shares[round][i][SHARES_TRIPLE+2], ((uint128_t *) MESSAGE1_DELTA_TRIPLE(message1) ) + round);
			}
		}

		// compute Legendre symbols of output
		for (int i = 0; i < RESIDUOSITY_SYMBOLS_PER_ROUND; ++i)
		{
			#ifdef LEGENDRE
				buffer1[ROUNDS*PARTIES*HASH_BYTES + round*RESIDUOSITY_SYMBOLS_PER_ROUND + i] = legendre_symbol_ct((uint128_t *) MESSAGE2_OUTPUT(message2) + round*RESIDUOSITY_SYMBOLS_PER_ROUND + i) ^ getBit(pk,((uint32_t *) challenge1)[round*RESIDUOSITY_SYMBOLS_PER_ROUND + i]);
			#else
				uint16_t prs = power_residue_symbol((uint128_t *) MESSAGE2_OUTPUT(message2) + round*RESIDUOSITY_SYMBOLS_PER_ROUND + i);
				prs += ((uint16_t) 254) - (uint16_t) pk[((uint32_t *) challenge1)[round*RESIDUOSITY_SYMBOLS_PER_ROUND + i]];
				prs %= 254;
				buffer1[ROUNDS*PARTIES*HASH_BYTES + round*RESIDUOSITY_SYMBOLS_PER_ROUND + i] = prs;
			#endif
		}	
	}

	unsigned char hash1[HASH_BYTES];
	HASH(buffer1,sizeof(buffer1),hash1);

	if( memcmp(hash1,MESSAGE1_COMMITMENT(message1),HASH_BYTES) != 0){
		printf("first hash fails!\n");
		return 0;
	}

	// check second commitment: alpha, beta and v
	uint128_t openings[ROUNDS][PARTIES][3] = {0};

	for (int round = 0; round < ROUNDS; ++round)
	{
		// import epsilon, alpha, beta and lambda
		uint128_t *epsilon = (uint128_t *) CHALLENGE2_EPSILON(challenge2) + round;
		uint128_t *alpha = (uint128_t *) MESSAGE3_ALPHA(message3)  + round;
		uint128_t *beta = (uint128_t *) MESSAGE3_BETA(message3)  + round;
		
		uint128_t sum_v_shares = 0;
		uint128_t sum_alpha_shares = 0;
		uint128_t sum_beta_shares = 0;
		uint128_t *lambda = (uint128_t *) CHALLENGE2_LAMBDA(challenge2) + round*RESIDUOSITY_SYMBOLS_PER_ROUND;

		// computes shares of alpha, beta and v
		for (int i = 0; i < PARTIES; ++i)
		{
			// TODO reuse code from respond 2
			if(i == unopened_party[round]){
				continue;
			}

			// compute share of alpha
			mul_add_mod_p(&openings[round][i][0],epsilon,&state.shares[round][i][SHARE_K]);
			add_mod_p(&openings[round][i][0],&state.shares[round][i][SHARES_TRIPLE]);
			reduce_mod_p(&openings[round][i][0]);
			add_mod_p(&sum_alpha_shares,&openings[round][i][0]);

			// compute share of beta and z
			uint128_t z_share = 0;
			openings[round][i][1] = state.shares[round][i][SHARES_TRIPLE+1];
			for (int j = 0; j < RESIDUOSITY_SYMBOLS_PER_ROUND; ++j)
			{
				// share of beta
				mul_add_mod_p(&openings[round][i][1],&state.shares[round][i][SHARES_R+j],lambda + j);

				// share of z
				uint128_t temp2 = 0;
				mul_add_mod_p(&temp2,&state.shares[round][i][SHARES_R + j],lambda + j); // this multiplication is done earlier, reuse result
				reduce_mod_p(&temp2);
				temp2 = m127 - temp2;
				uint128_t index = state.indices[round*RESIDUOSITY_SYMBOLS_PER_ROUND + j];
				mul_add_mod_p(&z_share,&temp2,&index); 

				if(i==0){
					mul_add_mod_p(&z_share, lambda+j , ((uint128_t *)MESSAGE2_OUTPUT(message2)) + round*RESIDUOSITY_SYMBOLS_PER_ROUND + j  );
				}
			}

			reduce_mod_p(&openings[round][i][1]);
			add_mod_p(&sum_beta_shares,&openings[round][i][1]);

			// compute sharing of v
			openings[round][i][2] = state.shares[round][i][SHARES_TRIPLE+ 2];
			if(i==0){
				mul_add_mod_p(&openings[round][i][2],alpha,beta);
			}
			reduce_mod_p(&openings[round][i][2]);
			openings[round][i][2] = m127 - openings[round][i][2];
			mul_add_mod_p(&openings[round][i][2],alpha,&state.shares[round][i][SHARES_TRIPLE+1]);
			mul_add_mod_p(&openings[round][i][2],beta,&state.shares[round][i][SHARES_TRIPLE]);
			mul_add_mod_p(&openings[round][i][2],epsilon,&z_share);
			reduce_mod_p(&openings[round][i][2]);
			add_mod_p(&sum_v_shares,&openings[round][i][2]);
		}

		// fill in unopened shares
		reduce_mod_p(&sum_alpha_shares);
		reduce_mod_p(&sum_beta_shares);
		reduce_mod_p(&sum_v_shares);
		openings[round][unopened_party[round]][0]  = m127 - sum_alpha_shares;
		openings[round][unopened_party[round]][1]  = m127 - sum_beta_shares;
		openings[round][unopened_party[round]][2]  = m127 - sum_v_shares;
		add_mod_p(&openings[round][unopened_party[round]][0],alpha);
		add_mod_p(&openings[round][unopened_party[round]][1],beta);
        reduce_mod_p(&openings[round][unopened_party[round]][0]);
        reduce_mod_p(&openings[round][unopened_party[round]][1]);
        reduce_mod_p(&openings[round][unopened_party[round]][2]);
	}

	unsigned char hash2[HASH_BYTES];
	HASH((unsigned char *)openings,sizeof(openings),hash2);

	if( memcmp(hash2,MESSAGE3_HASH(message3),HASH_BYTES) != 0){
		printf("second hash fails!\n");
		return 0;
	}

	return 1;
}
