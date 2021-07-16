#include <stdio.h>
#include <stdint.h>
#include "sign.h"
#include <time.h>

#define TRIALS 100

static inline
uint64_t rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((uint64_t)hi << 32) | lo;
}
#define TIC printf("\n"); uint64_t cl = rdtsc();
#define TOC(A) printf("%s cycles = %lu \n",#A ,rdtsc() - cl); cl = rdtsc();

int main(){

	//test();
	//return;

	int executions_done,i;
	unsigned char pk[PK_BYTES] = {0};
	unsigned char sk[SK_BYTES] = {0};

	printf("pk bytes : %d\n", PK_BYTES );
	printf("sk bytes : %d\n", SK_BYTES );
	printf("sig bytes : %ld\n", (uint64_t) SIG_BYTES );

	unsigned char message[1];
	message[0] = 42;
	unsigned char sig[SIG_BYTES];
	uint64_t sig_len;

	uint64_t keygenTime = 0;
	uint64_t signTime = 0;
	uint64_t verifyTime = 0;
	uint64_t t;

    clock_t keygenStart = 0;
    clock_t keygenEnd = 0;
    clock_t signStart = 0;
    clock_t signEnd = 0;
    clock_t verifyStart = 0;
    clock_t verifyEnd = 0;

    double keygenTotal = 0;
    double signTotal = 0;
    double verifyTotal = 0;


	for(int i=0 ; i<TRIALS; i++){
		if(i == 0){
            keygenStart = clock();
			t = rdtsc();
			keygen(pk,sk);
			keygenTime += rdtsc()-t;
            keygenEnd = clock();
		}


		signStart = clock();
		t = rdtsc();
		sign(sk,pk,message,1,sig,&sig_len);
		signTime += rdtsc()-t;
		signEnd = clock();

		verifyStart = clock();
		t = rdtsc();
		int ver = verify(pk,message,1,sig);
		verifyTime += rdtsc()-t;
		verifyEnd = clock();

        keygenTotal += (double)(keygenEnd - keygenStart);
        signTotal += (double)(signEnd - signStart);
        verifyTotal += (double)(verifyEnd - verifyStart);

		if(ver <= 0){
			printf("Signature invalid! \n");
		}
	}

	printf("keygen cycles :       %lu \n", keygenTime );
	printf("signing cycles :      %lu \n", signTime/TRIALS );
	printf("verification cycles : %lu \n", verifyTime/TRIALS );

    printf("keygen time :         %f ms\n", 1000*keygenTotal/CLOCKS_PER_SEC);
    printf("signing time :        %f ms\n", 1000*signTotal/(((double) TRIALS)*CLOCKS_PER_SEC));
    printf("verification time :   %f ms\n", 1000*verifyTotal/(((double) TRIALS)*CLOCKS_PER_SEC));

    return 0;
}
