#include "random.h"

/*
pRNG based on http://www.cs.wm.edu/~va/software/park/park.html
*/
#define MODULUS    2147483647
#define MULTIPLIER 48271
#define DEFAULT    123456789

static long seed = DEFAULT;

double Random(void)
/*
----------------------------------------------------------------
Random returns a pseudo-random real number uniformly distributed
between 0.0 and 1.0.
----------------------------------------------------------------
*/
{
	const long Q = MODULUS / MULTIPLIER;
	const long R = MODULUS % MULTIPLIER;
	long t;

	t = MULTIPLIER * (seed % Q) - R * (seed / Q);
	if (t > 0)
		seed = t;
	else
		seed = t + MODULUS;
	return ((double) seed / MODULUS);
}

/*
End of the pRNG algorithm
*/

