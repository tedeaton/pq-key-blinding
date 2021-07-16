/*
Implementation by the Keccak Team, namely, Guido Bertoni, Joan Daemen,
Michaël Peeters, Gilles Van Assche and Ronny Van Keer,
hereby denoted as "the implementer".

For more information, feedback or questions, please refer to our website:
https://keccak.team/

To the extent possible under law, the implementer has waived all copyright
and related or neighboring rights to the source code in this file.
http://creativecommons.org/publicdomain/zero/1.0/
*/

#include <string.h>
#include "test_crypto_hash.h"

int main()
{
    const unsigned char *message = (const unsigned char *)
        "\x00\x01\x02\x03\x04\x05\x06\x07\x08\x09\x0A\x0B\x0C\x0D\x0E\x0F\x10";
    const unsigned char *hash = (const unsigned char *)
        "\x6b\xf7\x5f\xa2\x23\x91\x98\xdb\x47\x72\xe3\x64\x78\xf8\xe1\x9b\x0f\x37\x12\x05\xf6\xa9\xa9\x3a\x27\x3f\x51\xdf\x37\x12\x28\x88\xb4\xb7\xa3\xa2\xb5\x98\xed\x4b\xd8\xfc\xf4\xcd\x38\xe0\x3d\xd8\x64\x74\xe4\x8e\xed\x8d\xb6\x41\x8d\xfe\x39\xa3\xd0\x7b\x85\x67\xb3\x33\x12\x65\xfa\x31\x03\x66\x6d\x9f\x62\xcd\x88\x7b\xfa\x4c\xf0\x74\x06\x46\xec\x8b\x8b\xc9\x17\xab\x2f\x22\x59\xc4\xb7\xcf\x7e\xfe\xd4\x41\x22\x5e\xd8\xa9\xad\xb0\x18\xa1\x7b\xca\xf4\xc8\xcb\x21\x94\x13\xcd\xc7\x37\xa9\x77\xfe\x59\xa6\x03\xc8\xc0\x02\x72\x9c\xe6\x7d\x44\x7e\x6e\xf7\xd5\x1c\xb1\xfb\x7b\x6a\x28\xf0\xad\x5b\x36\x2b\x45\x57\x66\x0b\xfb\x06\x2c\xc0\x4f\xd1\xc0\xe7\xb5\x3d\xc5\xd3\xa0\xfa\x7c\x54";
    return test_crypto_hash(message, 17, hash, 168);
}
