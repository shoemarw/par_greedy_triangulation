/*
 * par_greedy_triangle.c
 * 
 *
 * CS 470 Research Project.
 * Parallel version (using MPI).
 *
 * Name(s): Randy, Eliza
 */
 
#include <sys/time.h> 
#include "greedy_triangle.h"

// Timer marcros written by Professor Lam, borrowed from PA2
#define START_TIMER(NAME) gettimeofday(&tv, NULL); \
    double NAME ## _time = tv.tv_sec+(tv.tv_usec/1000000.0);
#define STOP_TIMER(NAME) gettimeofday(&tv, NULL); \
    NAME ## _time = tv.tv_sec+(tv.tv_usec/1000000.0) - (NAME ## _time);
#define GET_TIMER(NAME) (NAME##_time)

// Toggle to generate image. It is recommended that this be commented out
// if more than 1000 points are in the point set. This is because the 
// image created will only be useful for a point set of less than 1000.
//#define IMAGE


point_t* my_points; // a processes' portion of the point set.

int main(int argc, char *argv[]) {
	
	// Set up MPI stuff
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	
	// Make sure we get the expected input.
	if (argc != 2) {
		printf("Usage %s <filename>, argv[0] \n", argv[0]);
		exit(EXIT_FAILURE);
	}
	
	
	if (my_rank == 0) {
		// Open the input file for reading. 
		char *fn = argv[1];
		FILE* fin = open_file(fn, "r");
		
		  //                         //
		 //  Read points from file  //
		//                         // 
		
		// The first line of a file must contain a number indicating
		// the number of points in the file. Read this value and use
		// it to allocate storage for the points.
		long num_points;
		fscanf(fin, "%ld\n", &num_points);
		point_t* points = (point_t*) allocate(num_points * sizeof(point_t));
		
		// Read in and store the point s.
		double x, y;     // The Cartesian coordinates of a point.
		long i = 0;    // Index for storing points.

		while (fscanf(fin, "%lf %lf\n", &x, &y) == 2) {
			// Put the values in a point struct and store.
			point_t *p = (point_t*) allocate(sizeof(point_t));
			p->x = x;
			p->y = y;
			
			// Make sure input file didn't make num_points too small.
			if (i >= num_points) 
			{
				error("ERROR: the number of lines exceeds expectation\n");
			}
			
			points[i] = *p;
			i++;
			free(p);
		}
		fclose(fin);
	}
	
	
	// Scatter all the points in points[] among the processes and store
	// locally in my_points.
//TODO	MPI_Scatter(
	
	  //                      //
	 //  Generate all lines  //
	//                      //
	// utility struct for timing calls
    struct timeval tv;
	START_TIMER(generate)
	STOP_TIMER(generate)
	
	
	
	
	// Clean up and exit
	free(points);
	free(lines);
	MPI_Finalize();
	return (EXIT_SUCCESS);
}
