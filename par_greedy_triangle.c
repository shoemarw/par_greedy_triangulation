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

const int  TAG = 1;
const int  ROOT = 0;


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


	MPI_Aint point_type_extent;
	MPI_Type_get_extent(MPI_point_t, NULL, &point_type_extent);
	MPI_Aint disps_l[3];	
	disps_l[0] = 0;
	disps_l[1] = disps_l[0] + point_type_extent;
	disps_l[2] = disps_l[1] + point_type_extent;

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

	int send_counts[nprocs]; 				// an array of how many points each process will recieve / how many root sends
	int points_to_recv;						// a single number from send_counts
	int displs_point_scatter[nprocs];		// the displacements for the scatterv, significant only to root

	if (my_rank==0) {
		// use interger division to determin the base amount for points each process will recieve 
		long base_point_count = num_points/(long)nprocs;

		// get the remainder to see how many leftover points there are
		int remainder = num_points%nprocs;



		// fill the array with the base number, then if there are remainders left add one to the 
		// count of how many points the process will recieve.
		for (int i = 0; i < nprocs; i++) {
			send_counts[i] = base_point_count;
			if (remainder) {
				send_counts[i] += 1;
				remainder--;
			} // end if
		} // end for

		// build displacement array
		displs_point_scatter[0] = 0;
		for (int i = 1; i < nprocs; i++) {
			displs_point_scatter[i] = displs_point_scatter[i-1] + send_counts[i-1];
		}
		// send each process how many points it should expect
		MPI_Scatter(send_counts, 1, MPI_INT, &points_to_recv, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
		my_points = (point_t*) allocate(points_to_recv*sizeof(point_t));

		// send each process its points
		MPI_Scatterv(points, send_counts, displs_point_scatter, MPI_point_t, my_points, points_to_recv,
                 MPI_point_t, ROOT, MPI_COMM_WORLD);
	}
	else {
		MPI_Scatter(send_counts, 1, MPI_INT, &points_to_recv, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
		MPI_Scatterv(points, send_counts, displs_point_scatter, MPI_point_t, my_points, points_to_recv,
                 MPI_point_t, ROOT, MPI_COMM_WORLD);
	}



  // - - - - - - - //
 //  Eliza start  //
// - - - - - - - //
	int recv_buff; 					// Used to check how many objects will be sent in next MPI_send 
	int num_lines;					// Number of lines to be calculated
	line_t* recv_lines = 0; 		// Used by root only
	int recv_lines_count[nprocs];	// Used by root only
	int displs[nprocs];				// Used by root only
	line_t* my_lines;				// Array of the process's calculated lines
	long num_of_lines;				// number of lines


	num_lines = ((num_points)*(num_points-1))/2;
	my_lines = (line_t*) allocate(num_lines * sizeof(line_t));
	// calculate lines 

	for (int iteration_square = 1; iteration_square < nprocs; iteration_square *= 2) {
		if (my_rank&iteration_square) {
			int num_cur_point = sizeof(my_points)/sizeof(point_t);
			int send_to = (my_rank-iteration_square+nprocs)%nprocs;
			
			// send the number of points the receiver should expect
			MPI_Send(&num_cur_point, 1, MPI_INT, send_to, TAG, MPI_COMM_WORLD);
			// send the points
			MPI_Send(my_points, num_cur_point, MPI_point_t, send_to, TAG, MPI_COMM_WORLD);
			
			num_of_lines = sizeof(my_lines)/sizeof(line_t);
			// send the number of lines the receiver should expect
			MPI_Gather(&num_of_lines, 1, MPI_LONG, recv_lines_count, 1, MPI_LONG, 0, MPI_COMM_WORLD);
			// send the lines
			MPI_Gatherv(&my_lines, num_of_lines, MPI_line_t, recv_lines, recv_lines_count, 
					    displs, MPI_line_t, 0, MPI_COMM_WORLD);
			break;
		}
		else {
			int recv_from = my_rank+iteration_square; 

			// receive number of points
			MPI_Recv(&recv_buff, 1, MPI_LONG, recv_from, MPI_ANY_TAG, MPI_COMM_WORLD, 
					 MPI_STATUS_IGNORE);

			// receive points into new_points
			point_t* new_points = (point_t*) allocate(recv_buff*sizeof(point_t));
			MPI_Recv(new_points, recv_buff, MPI_point_t, recv_from, MPI_ANY_TAG, 
					 MPI_COMM_WORLD, MPI_STATUS_IGNORE);
		
			// calc new lines //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  


//resize my_points and add new_points	
		}
	}// end for

	if (my_rank==0) {
		// get the number of line_t each process is sending
		MPI_Gather(&num_of_lines, 1, MPI_INT, recv_lines_count, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

		displs[0] = 0;
		long total_line_num = recv_lines_count[0];

		// calculate how many total lines are being sent and the displs
        for (int i=1; i < nprocs; i++) {
           total_line_num += recv_lines_count[i];
           displs[i] = displs[i-1] + recv_lines_count[i-1];
        }

		recv_lines = (line_t*) allocate( total_line_num* sizeof(line_t));
		MPI_Gatherv(my_lines, num_of_lines, MPI_line_t, recv_lines, recv_lines_count, 
	    displs, MPI_line_t, ROOT, MPI_COMM_WORLD);


	    // root scatterv
	}
	else {
		// non-root scatterv
	}
  // - - - - - - - //
 //   Eliza end   //
// - - - - - - - //


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
	MPI_Type_free(&MPI_line_t);
	MPI_Finalize();
	return (EXIT_SUCCESS);
}
