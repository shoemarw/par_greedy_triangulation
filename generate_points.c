/*
 * generate_points.c
 * 
 * Creates a number of points in the standard x-y plane in
 * Cartesian coordinates. The number of points is supplied
 * via the command line as the first argument. The second
 * command line argument is the seed used to generate the 
 * pseudo-random numbers. A linear congruential generator 
 * is used to create the random numbers.
 *
 * Name(s): Randy and Eliza
 */

#include "greedy_triangle.h"

int main(int argc, char *argv[]) {
	
	// make sure we get the expected input.
	if (argc != 3) {
		printf("Usage %s <filename>, argv[0] <number of points>	<seed>\n",
				argv[0]);
		exit(EXIT_FAILURE);
	}
	
	// the number of points to be generated
	unsigned long num_points = strtol(argv[1], NULL, 10);
	// the seed for generating pseudo-random numbers
	long seed = strtol(argv[2], NULL, 10);
	
	// create the name of the file where the points will be stored.
	// the name will be "test<number of points>pts<seed>"
	char *fname;
	fname = allocate(strlen(argv[1]) + strlen(argv[2]) + 8);			 
	fname[0] = '\0';
	strcat(fname, "test");
	strcat(fname, argv[1]);
	strcat(fname, "pts");
	strcat(fname, argv[2]);
	
	// Create a file for writing
	FILE* fin = open_file(fname, "w");
	long long iterations = 2*num_points;
	
	// generate 2*num_points pseudo-random numbers
	char *pts = allocate(iterations);
	for (int i = 0; i < iterations; i++) {
		seed = (11*seed + 67) % 100;
		pts[i] = seed - 49;
	}
	
	// store the number of points on the first line of
	// the file
	fprintf(fin, "%s", argv[1]);
	fprintf(fin, "\n");
	// store the numbers in a file
	for (int i = 0; i < iterations; i += 2) {
		fprintf(fin, "%d", pts[i]);
		fprintf(fin, " ");
		fprintf(fin, "%d", pts[i+1]);
		fprintf(fin, "\n");
	}
	
	fclose(fin);
	free(fname);
	free(pts);
	return (EXIT_SUCCESS);
}