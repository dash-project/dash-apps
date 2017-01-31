#ifndef PACKING_DNA_SEQ_H
#define PACKING_DNA_SEQ_H

#include <math.h>
#include <assert.h>
#include <string.h>

unsigned int packedCodeToFourMer[256];

#define pow4(a) (1<<((a)<<1))

void init_LookupTable()
{
	// Work with 4-mers for the moment to have small lookup tables
	int merLen = 4, i, slot, valInSlot;
	unsigned char mer[4];

	for ( i = 0; i < 256; i++ ) {
		// convert a packedcode to a 4-mer
		int remainder = i;
		int pos = 0;
		for( slot = merLen-1; slot >= 0; slot-- ) {
         valInSlot = remainder / pow4(slot);
			char base;

			if( valInSlot == 0 ) { base = 'A'; }
			else if( valInSlot == 1 ) { base = 'C'; }
			else if( valInSlot == 2 ) { base = 'G'; }
			else if( valInSlot == 3 ) { base = 'T'; }
			else{ assert( 0 ); }

			mer[pos] = base;
			pos++;
         remainder -= valInSlot * pow4(slot);
		}
		unsigned int *merAsUInt = (unsigned int*) mer;
		packedCodeToFourMer[i] = (unsigned int) (*merAsUInt);
	}
}

unsigned char convertFourMerToPackedCode(unsigned char *fourMer)
{
	int retval = 0;
	int code, i;
	int pow = 64;

	for ( i=0; i < 4; i++) {
		char base = fourMer[i];
		switch ( base ) {
			case 'A':
				code = 0;
				break;
			case 'C':
				code = 1;
				break;
			case 'G':
				code = 2;
				break;
			case 'T':
				code = 3;
				break;
		}
		retval += code * pow;
		pow /= 4;
	}
	return ((unsigned char) retval);
}

void packSequence(const unsigned char *seq_to_pack, unsigned char *m_data, int m_len)
{
	/* The pointer to m_data points to the result of the packing */

	int ind, j = 0;     // coordinate along unpacked string ( matches with m_len )
    int i = 0;     // coordinate along packed string

    // do the leading seq in blocks of 4
    for ( ; j <= m_len - 4; i++, j+=4 ) {
		m_data[i] = convertFourMerToPackedCode( ( unsigned char * ) ( seq_to_pack + j )) ;
    }

    // last block is special case if m_len % 4 != 0: append "A"s as filler
    int remainder = m_len % 4;
    unsigned char blockSeq[5] = "AAAA";
    for(ind = 0; ind < remainder; ind++) {
		blockSeq[ind] = seq_to_pack[j + ind];
    }
    m_data[i] = convertFourMerToPackedCode(blockSeq);
}

void unpackSequence(const unsigned char *seq_to_unpack, unsigned char *unpacked_seq, int kmer_len)
{
	/* Result string is pointer unpacked_seq */
	int i = 0, j = 0;
   int packed_len = (kmer_len+3)/4;
	for( ; i < packed_len ; i++, j += 4 )
	{
		*( ( unsigned int * )( unpacked_seq + j ) ) = packedCodeToFourMer[ seq_to_unpack[i] ];
	}
	*(unpacked_seq + kmer_len) = '\0';
}

int comparePackedSeq(const unsigned char *seq1, const unsigned char *seq2, int seq_len)
{
	return memcmp(seq1, seq2, seq_len);
}

#endif // PACKING_DNA_SEQ_H

