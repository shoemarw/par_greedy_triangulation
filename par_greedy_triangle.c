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
#include <mpi.h>

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

int nprocs;         // number of processes
int my_rank;        // rank of a particular process
point_t* my_points; // a processes' portion of the point set.
line_t*  my_lines;  // a processes' portion of lines.

//TEMPORARILY MAKING THESE GLOBAL, THIS MAY NEED TO BE CHANGED
long num_points;
point_t* points;
line_t* lines;


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

	// Create MPI derived data types needed to communicate points and lines.
	// We time this to see how much overhead cost it adds in terms of time.
	// utility struct for timing calls
    struct timeval tv;
	START_TIMER(MPIoverhead)

	//Build MPI_point_t
	MPI_Datatype MPI_point_t, MPI_line_t;

	int block_lens_p[] = {2};

	MPI_Aint disps_p[2];
	disps_p[1] = offsetof(point_t, x);
	disps_p[2] = offsetof(point_t, y);

	MPI_Datatype types_p[2] = {MPI_DOUBLE,MPI_DOUBLE};

	MPI_Type_create_struct(2, block_lens_p, disps_p, types_p, &MPI_point_t);
	MPI_Type_commit(&MPI_point_t);

	//Build MPI_line_t
	int block_lens_l[] = {1,1,1};

	MPI_Aint disps_l[3];	
	disps_l[0] = offsetof(line_t, p);
	disps_l[1] = offsetof(line_t, q);
	disps_l[2] = offsetof(line_t, len);

	MPI_Datatype types_l[3] = {MPI_point_t, MPI_point_t, MPI_DOUBLE};
	
	MPI_Type_create_struct(3, block_lens_l, disps_l, types_l, &MPI_line_t);

	MPI_Type_commit(&MPI_line_t);

	STOP_TIMER(MPIoverhead)

	

	
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
	
	  //                      //
	 //  Generate all lines  //
	//                      //
	
	START_TIMER(generate)
	// Scatter all the points in points[] among the processes and store
	// locally in my_points so Eliza's binomial tree structure can be used.
//TODO	MPI_Scatter(

//////////////////////////////////////////////////////////////////////////////
//TEMPORARY SERIAL CODE SO THAT PHASES 2 AND 3 CAN BE DEVELOPED...////////////
//THIS SHOULD BE REMOVED WHEN ELIZA'S BINOMIAL TREE STRUCTURE IS IMPLEMENTED//
	if (my_rank == 0) {														//
		int num_lines = ((num_points)*(num_points-1))/2;					//
		line_t* lines = (line_t*) allocate(num_lines * sizeof(line_t));		//
																			//
		long index = 0;														//
		for (int i = 0; i < num_points; i++) {								//
			// Compute the distance between point i and every point			//
			// from i+1 onward. Then 'make' and store the corresponding		//
			// line.														//
			for (int j = i+1; j < num_points; j++) {						//
				double length = distance(&points[i], &points[j]);			//
				line_t* l = (line_t*) allocate(sizeof(line_t));				//
				// set the values of the line and store it.					//
				l->p =         &points[i];									//
				l->q =         &points[j];									//
				l->len =       length;										//
				lines[index] = *l;											//
				index++;													//
				free(l);													//
			}																//
		}																	//
	}																		//
//TEMPORARY CODE TO SCATTER THE LINES FROM PROC 0 TO ALL PROCS////////////////
//	MPI_Scatter(lines);
//////////////////////////////////////////////////////////////////////////////


	STOP_TIMER(generate)
	
	  //                                      //
	 //  Sort the lines from small to large  //
	//                                      //
	
	START_TIMER(sort)
//	qsort(lines, num_lines, sizeof(line_t), compare);
	STOP_TIMER(sort)
	
	  //                                   //
     //  Greedily build the tringulation  //
    //	                                 //
	
	START_TIMER(triangulate)
	STOP_TIMER(triangulate)
	
	  //                                     //
     //  Triangulation Done, Display Stats  //
    //	                                   //
	
	// These stats are only for the portions of the code specific to the three
	// phases of building the greedy triangulation. Generate all lines, sort the
	// lines in non-decreasing order, and greedily adding line segments to the
	// triangulation.
	printf("Gent: %.4f  Sort: %.4f  Tria: %.4f\n Overhead: %.4f\n",
	        GET_TIMER(generate), GET_TIMER(sort), 
	        GET_TIMER(triangulate), GET_TIMER(MPIoverhead));
	
	
	// Clean up and exit
	free(points);
	free(lines);
	MPI_Type_free(&MPI_point_t);
	// MPI_Type_free(&MPI_line_t);
	MPI_Finalize();
	return (EXIT_SUCCESS);
}
