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

	uint8_t *marked = new uint8_t[size >> 1];
	if (marked == NULL) {
		printf("Cannot allocate enough memory\n");
		MPI_Finalize();
		return 1;
	}

	for (uint64_t i = 0; i < (size >> 1); i++)
		marked[i] = 0;
	//marked[0] was for 1, but since 1 is not prime but 2 is, so marked[0] is for 2 instead.

	//First Prime starts at 3
	uint64_t prime = 3;
	uint64_t first;
	while (prime * prime <= n) {
		if (prime * prime > low_value)
			first = prime * prime - low_value;
		else {
			if (low_value % prime == 0)
				first = 0;
			else
				first = prime - (low_value % prime);
		}
		if (((low_value + first) & 1) == 0) {
			first += prime;
		}
		for (uint64_t i = first; i < size; i += (prime + prime))
			if (i & 1)//only if i is odd
				marked[i >> 1] = 1;
		if (id == 0) {
			uint64_t next_Unmark = (prime >> 1) + 1;
			while ((marked[next_Unmark]) && (next_Unmark < (high_value >> 1))) {
				next_Unmark++;
			}
			//next prime
			prime = (next_Unmark << 1) + 1; //*2+1, same speed if -O3.
		}
		MPI_Bcast(&prime, 1, MPI_UINT64_T, 0, MPI_COMM_WORLD);
	}


	int count = 0;
	for (uint64_t i = 0; i < (size >> 1); i++)
		if (!marked[i])
			count++;


	int global_count = 0;
	MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	elapsed_time += MPI_Wtime();


	if (id == 0) {
		printf("The total number of prime: %d, total time: %10.6f, total node ", global_count, elapsed_time);
		//printf("%d primes are less than or equal to %llu\n",
		//	global_count, n);
		//printf("Total elapsed time: %10.6f\n", elapsed_time);
		//printf("This is Part 1\n");
	}

	MPI_Finalize();
	return 0;
}
