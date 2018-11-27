#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>

#define MIN(a,b)  ((a)<(b)?(a):(b))
#define BLOCK_LOW(id,p,n)   ((id)*(n)/(p))
#define BLOCK_HIGH(id,p,n)  (BLOCK_LOW((id)+1,p,n)-1)
#define BLOCK_SIZE(id,p,n)  \
            (BLOCK_HIGH(id,p,n)-BLOCK_LOW(id,p,n)+1)


long long atoll(const char *p){
	long long n;
	int c, neg = 0;
	unsigned char   *up = (unsigned char *)p;

	if (!isdigit(c = *up)) {
		while (isspace(c))
			c = *++up;
		switch (c) {
		case '-':
			neg++;
		case '+':
			c = *++up;
		}
		if (!isdigit(c))
			return (0);
	}

	for (n = '0' - c; isdigit(c = *++up); ) {
		n *= 10; /* two steps to avoid unnecessary overflow */
		n += '0' - c; /* accum neg to avoid surprises at MAX */
	}

	return (neg ? n : -n);
}

int main (int argc, char *argv[])
{
	double elapsed_time;
	long long low_value, high_value, proc0_size;
	int id, p;
	char* marked;
	long long index;
	long long first,i;
	long long n,size;
	long long global_count, count, prime;
   long long node;
	MPI_Init (&argc, &argv);
 	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = -MPI_Wtime();
	MPI_Comm_rank (MPI_COMM_WORLD, &id);
   	MPI_Comm_size (MPI_COMM_WORLD, &p);
	if (argc != 3) {
     	if (!id) printf ("Command line: %s <m>\n", argv[0]);
      	MPI_Finalize(); 
      	exit (1);
	}
	n = atoll(argv[1]);
   node = atoll(argv[2]);
   	low_value = 2 + BLOCK_LOW(id,p,n-1);
   	high_value = 2 + BLOCK_HIGH(id,p,n-1);
   	size = BLOCK_SIZE(id,p,n-1);

      if(low_value %2 == 0) {
         if (high_value % 2 == 0) {
            size = (long long)floor(size / 2.0);
            -- high_value; 
         }
         else {
            size = size/2;
         }
         ++ low_value;
      }
      else {
         if(high_value % 2 == 0) {
            size = size/2;
            -- high_value;
         }
         else size = (long long)ceil(size/2.0);
      }

   	proc0_size = (n-1)/p;
   	if ((2 + proc0_size) < (long long) sqrt((double) n)) {
      	if (!id) printf ("Too many processes\n");
      	MPI_Finalize();
      	exit (1);
   	}

   	marked = (char *) malloc (size);
   	if (marked == NULL) {
      	printf ("Cannot allocate enough memory\n");
      	MPI_Finalize();
      	exit (1);
   	}
   	for (i = 0; i < size; i++) marked[i] = 0;
   	if (!id) index = 0;
   	prime = 3;
   	do {
      	if (prime * prime > low_value)
         	first = (prime * prime - low_value)/2;
      	else {
         	if (!(low_value % prime)) first = 0;
         	else {
               if((low_value/prime) % 2 != 0) {
                  first = (2*prime - (low_value % prime)) / 2;
               }
               else {
                  first = (prime - (low_value % prime)) / 2;
               }
            }
      	}
      	for (i = first; i < size; i += prime) marked[i] = 1;
      	if (!id) {
         	while (marked[++index]);
         	prime = index * 2 + 3;
      	}
      	MPI_Bcast (&prime,  1, MPI_INT, 0, MPI_COMM_WORLD);
   	} while (prime * prime <= n);
   	count = 0;
   	for (i = 0; i < size; i++)
      	if (!marked[i]) count++;
   	MPI_Reduce (&count, &global_count, 1, MPI_INT, MPI_SUM,0, MPI_COMM_WORLD);
   	elapsed_time += MPI_Wtime();
   	if (!id) {
         ++ global_count;
      	printf ("The total number of prime: %lld, total time: %10.6f , total node %d\n",global_count, elapsed_time, node);
   	}
   	MPI_Finalize ();
   	return 0;
}
