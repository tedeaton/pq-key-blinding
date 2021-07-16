CC=gcc
CFLAGS=-I XKCP/bin/generic64/
LFLAGS=-L XKCP/bin/generic64/ -lkeccak -lcrypto -L/usr/local/opt/openssl/lib -I/usr/local/opt/openssl/include

IMPLEMENTATION_SOURCE = sign.c test.c zkproof.c merkletree.c
IMPLEMENTATION_HEADERS= sign.h keccaklib zkproof.h parameters.h merkletree.h

BLIND_SOURCE = blsign.c bltest.c blzkproof.c merkletree.c
BLIND_HEADERS = blsign.h keccaklib blzkproof.h parameters.h merkletree.h

test: $(IMPLEMENTATION_SOURCE) $(IMPLEMENTATION_HEADERS)
	gcc -o test $(IMPLEMENTATION_SOURCE) $(CFLAGS) $(LFLAGS) -std=c11 -O3 -g -march=native 

blindtest: $(BLIND_SOURCE) $(BLIND_HEADERS) 
	gcc -o blindtest $(BLIND_SOURCE) $(CFLAGS) $(LFLAGS) -std=c11 -O3 -g -march=native

keccaklib: 
	(cd XKCP; make generic64/libkeccak.a)

.PHONY: clean
clean:
	rm -f PQCgenKAT_sign test debug test_offline intermediateValues.txt *.req *.rsp >/dev/null
