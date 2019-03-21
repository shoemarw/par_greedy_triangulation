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



void read_points(char *argv[]) {

	// Open the input file for reading. 
	char *fn = argv[1];
	FILE* fin = open_file(fn, "r");
	

	// The first line of a file must contain a number indicating
	// the number of points in the file. Read this value and use
	// it to allocate storage for the points.
	fscanf(fin, "%ld\n", &num_points);
	points = (point_t*) allocate(num_points * sizeof(point_t));
	
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


void distrib_points() {

	int send_counts[nprocs];			// an array of how many points each process will recieve / how many root sends
	int displs_point_scatter[nprocs];	// the displacements for the scatterv, significant only to root
	int points_to_recv = 0;				// used to store a single number from send_counts

	if (my_rank==ROOT) {
		// use interger division to determine the base amount for points each process will recieve 
		long base_point_count = num_points/(long)nprocs;

		// get the remainder to see how many leftover points there are
		int remainder = num_points%nprocs;


		// fill the array with the base number, then if there are remainders left add one to the 
		// count of how many points the process will recieve.
		for (int i = 0; i < nprocs; i++) {
			send_counts[i] = base_point_count * sizeof(point_t);
			if (remainder > 0) {
				send_counts[i] += 1 * sizeof(point_t);
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


		my_points = (point_t*) allocate((int)(points_to_recv*sizeof(point_t)));
		// send each process its points
		MPI_Scatterv(points, send_counts, displs_point_scatter, MPI_BYTE, my_points, points_to_recv,
                 MPI_BYTE, ROOT, MPI_COMM_WORLD);

		printf("%lf %lf\n", my_points[0].x, my_points[0].y);
	}
	else { // NOT root
		MPI_Scatter(send_counts, 1, MPI_INT, &points_to_recv, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

		my_points = (point_t*) allocate((int)(points_to_recv*sizeof(point_t)));
		MPI_Scatterv(points, send_counts, displs_point_scatter, MPI_BYTE, my_points, points_to_recv,
                 MPI_BYTE, ROOT, MPI_COMM_WORLD);
	}
}



void gen_lines() {

	// num lines = points*(points-1)/2
	num_points = sizeof(my_points)/sizeof(point_t);
	int num_lines = ((num_points)*(num_points-1))/2;
	my_lines = (line_t*) allocate(num_lines * sizeof(line_t));
	
	long index = 0;
	for (int i = 0; i < num_points; i++) {
		// Compute the distance between point i and every point
		// from i+1 onward. Then 'make' and store the corresponding
		// line.
		for (int j = i+1; j < num_points; j++) {
			double length = distance(&points[i], &points[j]);
			line_t* l = (line_t*) allocate(sizeof(line_t));
			// set the values of the line and store it.
			l->p =         &points[i];
			l->q =         &points[j];
			l->len =       length;
			my_lines[index] = *l;
			index++;
			free(l);
		}
	}

	int recv_buff; 	// Used to check how many objects will be sent in MPI_send 
	// In this for loop we calculated all the lines
	for (int iteration_square = 1; iteration_square < nprocs; iteration_square *= 2) {
		if (my_rank&iteration_square) {
			long sizeof_cur_points = sizeof(my_points);
			int send_to = (my_rank-iteration_square+nprocs)%nprocs;

			// send the number of points the receiver should expect
			MPI_Send(&sizeof_cur_points, 1, MPI_LONG, send_to, TAG, MPI_COMM_WORLD);
			// send the points
			MPI_Send(my_points, sizeof_cur_points, MPI_BYTE, send_to, TAG, MPI_COMM_WORLD);
			
			break;
		}
		else {
			int recv_from = my_rank+iteration_square; 
			// receive number of points
			MPI_Recv(&recv_buff, 1, MPI_LONG, recv_from, MPI_ANY_TAG, MPI_COMM_WORLD, 
					 MPI_STATUS_IGNORE);

			// receive points into new_points
			point_t* new_points = (point_t*) allocate(recv_buff);
			MPI_Recv(new_points, recv_buff, MPI_BYTE, recv_from, MPI_ANY_TAG, 
					 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			int num_my_point = sizeof(my_points)/sizeof(point_t);
			int num_new_points = sizeof(new_points)/sizeof(point_t);
			line_t* new_lines = (line_t*) allocate(sizeof(line_t)*num_my_point*num_new_points);
			int new_line_index = 0;
			for (int j = 0; j < num_my_point; j++){
				for (int k = 0; k < num_new_points; ++k)
				{
					double length = distance(&my_points[j], &new_points[k]);
					line_t* l = (line_t*) allocate(sizeof(line_t));
					l->p =         &my_points[j];
					l->q =         &new_points[k];
					l->len =       length;
					new_lines[new_line_index] = *l;
					new_line_index++;
					free(l);
				}
			}
			line_t* temp = (line_t*) allocate(sizeof(my_lines)+sizeof(new_lines));
			int num_my_linr = sizeof(my_lines)/sizeof(line_t);
			copy_array(my_lines, temp, num_my_linr);
			copy_array(new_lines, temp[num_my_linr], sizeof(new_lines)/sizeof(line_t))

		}
	}// end for
}


void distrib_lines() {
	long num_of_lines;				// number of lines
	int displs[nprocs];				// Used by root only
	line_t* recv_lines = 0; 		// Used by root only
	int* recv_lines_count;			// Used by root only

	if (my_rank == 0) {
		recv_lines_count = allocate (sizeof(int) * nprocs);
	}

	num_of_lines = sizeof(my_lines);
	// send the number of lines the receiver should expect
	MPI_Gather(&num_of_lines, 1, MPI_INT, recv_lines_count, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

	if (my_rank==ROOT) {
		displs[0] = 0;
		long total_line_num = recv_lines_count[0];

		// calculate how many total lines are being sent and the displs
        for (int i=1; i < nprocs; i++) {
           total_line_num += recv_lines_count[i];
           displs[i] = displs[i-1] + recv_lines_count[i-1];
        }
		recv_lines = (line_t*) allocate( total_line_num* sizeof(line_t));        
	}

	// send all lines to ROOT
	MPI_Gatherv(&my_lines, num_of_lines, MPI_BYTE, recv_lines, recv_lines_count, 
			    displs, MPI_BYTE, ROOT, MPI_COMM_WORLD);

//  //  //  //  //  // scatter lines //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  //  
}




////			  //
/// 	Main 	 ///
// 				////

int main(int argc, char *argv[]) {

	// Set up MPI stuff
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);	
MPI_Barrier(MPI_COMM_WORLD); if (my_rank==ROOT) printf("line 273, nprocs %i\n", nprocs);
	
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

	STOP_TIMER(MPIoverhead)

	

MPI_Barrier(MPI_COMM_WORLD); if (my_rank==ROOT) printf("line 267\n");
	// Root reads in the lines for given file
	if (my_rank==ROOT) {
		read_points(argv);
	}

MPI_Barrier(MPI_COMM_WORLD); if (my_rank==ROOT) printf("line 273\n");

	// Root scatters the points
	distrib_points();
MPI_Barrier(MPI_COMM_WORLD); if (my_rank==ROOT) printf("line 277\n");


	
	START_TIMER(generate)

	// gen_lines();

	// distrib_lines();

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
	if (my_rank==ROOT) {
		printf("Gent: %.4f  Sort: %.4f  Tria: %.4f\n Overhead: %.4f\n",
	        GET_TIMER(generate), GET_TIMER(sort), 
	        GET_TIMER(triangulate), GET_TIMER(MPIoverhead));
	}
	
	// Clean up and exit
	// free(points);
	// free(lines);
	MPI_Finalize();
	return (EXIT_SUCCESS);
}
