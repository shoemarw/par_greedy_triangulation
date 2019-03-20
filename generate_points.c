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
	
	// Make sure we get the expected input.
	if (argc != 3) {
		printf("Usage %s <filename>, argv[0] <number of points>	<seed>\n",
				argv[0]);
		exit(EXIT_FAILURE);
	}
	
	// The number of points to be generated
	unsigned long num_points = strtol(argv[1], NULL, 10);
	// The seed for generating pseudo-random numbers
	unsigned long seed = strtol(argv[2], NULL, 10);
	
	// Set the range in which all x-y values are generated
	int range;
	if (num_points < 10000) {
		range = 100;
	} else {
		range = sqrt(num_points);
	}
	
	// Create the name of the file where the points will be stored.
	// The name will be "test<number of points>pts<seed>"
	char *fname;
	fname = allocate(strlen(argv[1]) + strlen(argv[2]) + 8);			 
	fname[0] = '\0';
	strcat(fname, "test");
	strcat(fname, argv[1]);
	strcat(fname, "pts");
	strcat(fname, argv[2]);
	
	// Create a file for writing
	FILE* fin = open_file(fname, "w");
		
	// Generate 2*num_points pseudo-random numbers
	char *pts = allocate(2*num_points);
	for (int i = 0; i < 2*num_points; i++) {
		// Compute the next pseudo random number in the sequence
		seed = (1103515245*seed + 12345) % ((1<<31));
		// Mod the value into the proper range and shift it so the
		// values fall into the interval [-(range/2)-1, range/2]
		pts[i] = seed%range - (range/2 +1);
	}
	
	// Store the number of points on the first line of the file
	fprintf(fin, "%s", argv[1]);
	fprintf(fin, "\n");
	// store the numbers in the file
	for (int i = 0; i < 2*num_points; i += 2) {
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