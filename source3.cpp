#include "mpi.h"
#include <cmath>
#include <cstdio>
#include <cstring>
#include <cstdint>
#include <cstdlib>
#define MIN(a,b)  ((a)<(b)?(a):(b))

#define BLOCK_LOW(id,p,n)  ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n) \
        ( BLOCK_LOW((id)+1,p,n)-1 ) 
#define BLOCK_SIZE(id,p,n) \
        (BLOCK_LOW( (id)+1, p, n) -  BLOCK_LOW((id), p, n))

#define BLOCK_OWNER(index,p,n) \
        ( ( ((p)*(index)+1)-1 ) / (n) )

void myBitSet(uint8_t *a, uint64_t pos) {
	a[pos >> 3] |= 1 << (pos & 0x07);
}

bool myBitCheck(uint8_t *a, uint64_t pos) {
	return a[pos >> 3] & (1 << (pos & 0x07));
}

int main(int argc, char *argv[])
{
	int rc = MPI_Init(&argc, &argv);
	MPI_Barrier(MPI_COMM_WORLD);
	double elapsed_time = -MPI_Wtime();
	int id = 0, p = 1;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	if (argc != 2) {
		if (id == 0)
			printf("Command line: %s <m>\n", argv[0]);
		MPI_Finalize();
		return 1;
	}
	//NO OFFSET
	uint64_t n = atoll(argv[1]);
	uint64_t low_value = BLOCK_LOW(id, p, n);
	uint64_t high_value = BLOCK_HIGH(id, p, n);
	uint64_t size = BLOCK_SIZE(id, p, n);
	uint64_t proc0_size = n / p;
	if (proc0_size < sqrt(n)) {
		if (id == 0)
			printf("Too many processes\n");
		MPI_Finalize();
		return 1;
	}

	uint8_t *marked = new uint8_t[(size >> 4) + 1];
	uint64_t numOfSP = sqrt(n);
	numOfSP = numOfSP >> 1;
	uint8_t *my_SP = new uint8_t[(numOfSP >> 3) + 1];
	if ((marked == NULL) || (my_SP == NULL)) {
		printf("Cannot allocate enough memory\n");
		MPI_Finalize();
		return 1;
	}
	for (uint64_t i = 0; i < (size >> 4) + 1; i++)
		marked[i] = 0;
	//marked[0] was for 1, but since 1 is not prime but 2 is, so marked[0] is for 2 instead.
	for (uint64_t i = 0; i < (numOfSP >> 3) + 1; i++)
		my_SP[i] = 0;


	uint64_t maxN = sqrt(numOfSP * 2);
	for (uint64_t i = 3; i <= maxN; i += 2) {
		if (myBitCheck(my_SP, i >> 1) == 0) {
			uint64_t j = i + i + i;
			while (j <= numOfSP) {
				myBitSet(my_SP, j >> 1);
				j += i + i;
			}
		}
	}



	//First Prime starts at 3
	uint64_t prime;
	uint64_t first;
	const uint64_t BlockSize = 1ULL << 21;

	for (uint64_t p = 0; p < size; p += BlockSize) {
		prime = 3;
		uint64_t blockLow = low_value + p;
		uint64_t blockHigh = MIN(blockLow + BlockSize, high_value);
		uint64_t realBlockSize = blockHigh - blockLow + 1;
		while (prime * prime <= blockHigh) {
			if (prime * prime > blockLow)
				first = prime * prime - blockLow;
			else {
				if (blockLow % prime == 0)
					first = 0;
				else
					first = prime - (blockLow % prime);
			}
			if (((blockLow + first) & 1) == 0) {
				first += prime;
			}
			for (uint64_t i = first; i < realBlockSize; i += (prime + prime))
				myBitSet(marked, (i+p) >> 1);
			uint64_t next_Unmark = (prime >> 1) + 1;
			while (myBitCheck(my_SP, next_Unmark) && (next_Unmark < numOfSP)) {
				next_Unmark++;
			}
			//next prime
			if (next_Unmark == numOfSP) {
				break;
			}
			else {
				prime = (next_Unmark << 1) + 1; //*2+1, same speed if -O3.
			}
		}
	}


	int count = 0;
	for (uint64_t i = 0; i < (size >> 1); i++)
		if (!myBitCheck(marked,i))
			count++;


	int global_count = 0;
	MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	elapsed_time += MPI_Wtime();


	if (id == 0) {
		printf("The total number of prime: %d, total time: %10.6f, total node ", global_count, elapsed_time);
		//printf("%d primes are less than or equal to %llu\n",
		//	global_count, n);
		//printf("Total elapsed time: %10.6f\n", elapsed_time);
		//printf("This is Part 3\n");
	}

	MPI_Finalize();
	return 0;
}
