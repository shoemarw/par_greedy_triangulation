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
point_t* pt_my_points; 	// a processes' portion of the point set.
line_t*  ln_my_lines;  	// a processes' portion of lines.
double* d_my_lines;		// a processes' portion of lines while in double format.

//TEMPORARILY MAKING THESE GLOBAL, THIS MAY NEED TO BE CHANGED
long l_num_points;
point_t* points;
line_t* lines;



void double_array_to_struct(double* arr, long size){
	//do stuff
}

void read_points(char *argv[]) {

	// Open the input file for reading. 
	char *fn = argv[1];
	FILE* fin = open_file(fn, "r");
	

	// The first line of a file must contain a number indicating
	// the number of points in the file. Read this value and use
	// it to allocate storage for the points.
	fscanf(fin, "%ld\n", &l_num_points);
	points = (point_t*) allocate(l_num_points * sizeof(point_t));
	
	// Read in and store the point s.
	double x, y;     // The Cartesian coordinates of a point.
	long i = 0;    // Index for storing points.

	while (fscanf(fin, "%lf %lf\n", &x, &y) == 2) {
		// Put the values in a point struct and store.
		point_t *p = (point_t*) allocate(sizeof(point_t));
		p->x = x;
		p->y = y;
		
		// Make sure input file didn't make l_num_points too small.
		if (i >= l_num_points) 
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

	int i_send_counts[nprocs];	// an array of how many points each process will recieve / how many root sends
	int i_displs_p[nprocs];		// the displacements for the scatterv, significant only to root
	int i_points_to_recv = 0;	// used to store a single number from i_send_counts

	if (my_rank==ROOT) {
		// use interger division to determine the base amount for points each process will recieve 
		long base_point_count = l_num_points/(long)nprocs;

		// get the remainder to see how many leftover points there are
		int remainder = l_num_points%nprocs;


		// fill the array with the base number, then if there are remainders left add one to the 
		// count of how many points the process will recieve.
		for (int i = 0; i < nprocs; i++) {
			i_send_counts[i] = base_point_count * sizeof(point_t);
			if (remainder > 0) {
				i_send_counts[i] += 1 * sizeof(point_t);
				remainder--;
			} // end if
		} // end for

		// build displacement array
		i_displs_p[0] = 0;
		for (int i = 1; i < nprocs; i++) {
			i_displs_p[i] = i_displs_p[i-1] + i_send_counts[i-1];
		}		
		// send each process how many points it should expect
		MPI_Scatter(i_send_counts, 1, MPI_INT, &i_points_to_recv, 1, MPI_INT, ROOT, MPI_COMM_WORLD);


		pt_my_points = (point_t*) allocate((int)(i_points_to_recv*sizeof(point_t)));
		// send each process its points
		MPI_Scatterv(points, i_send_counts, i_displs_p, MPI_BYTE, pt_my_points, i_points_to_recv,
                 MPI_BYTE, ROOT, MPI_COMM_WORLD);
	}
	else { // NOT root
		MPI_Scatter(i_send_counts, 1, MPI_INT, &i_points_to_recv, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

		pt_my_points = (point_t*) allocate((int)(i_points_to_recv*sizeof(point_t)));
		MPI_Scatterv(points, i_send_counts, i_displs_p, MPI_BYTE, pt_my_points, i_points_to_recv,
                 MPI_BYTE, ROOT, MPI_COMM_WORLD);
	}
}



void gen_lines() {

	// Generate lines as array of (5) doubles		
	l_num_points = sizeof(pt_my_points)/sizeof(point_t);
	int i_num_lines = ((l_num_points)*(l_num_points-1))/2;
	d_my_lines = (double*) allocate(i_num_lines * sizeof(double) * 5);
	const int X0 = 0;
	const int Y0 = 1;
	const int X1 = 2;
	const int Y1 = 3;
	const int LEN = 4;
	long l_index = 0;


	for (int i = 0; i < l_num_points; i++) {
		// Compute the distance between point i and every point
		// from i+1 onward. Then 'make' and store the corresponding
		// line.
		for (int j = i+1; j < l_num_points; j++) {
			double length = distance(&pt_my_points[i], &pt_my_points[j]);
			d_my_lines[l_index + X0] = pt_my_points[i].x;
			d_my_lines[l_index + Y0] = pt_my_points[i].y;
			d_my_lines[l_index + X1] = pt_my_points[j].x;
			d_my_lines[l_index + Y1] = pt_my_points[j].y;
			d_my_lines[l_index + LEN] = length;

			l_index +=5;
		}
	}


	int recv_buff; 	// Used to check how many objects will be sent in MPI_send 
	// In this for loop we calculated all the lines
	for (int iteration_square = 1; iteration_square < nprocs; iteration_square *= 2) {
		if (my_rank&iteration_square) {
			long l_size_cur_pt = sizeof(pt_my_points);
			int i_send_to = (my_rank-iteration_square+nprocs)%nprocs;

			// send the number of points the receiver should expect
			MPI_Send(&l_size_cur_pt, 1, MPI_LONG, i_send_to, TAG, MPI_COMM_WORLD);
			// send the points
			MPI_Send(pt_my_points, l_size_cur_pt, MPI_BYTE, i_send_to, TAG, MPI_COMM_WORLD);
			
			break;
		}
		else {
			int i_recv_from = my_rank+iteration_square; 
			// receive number of points
			MPI_Recv(&recv_buff, 1, MPI_LONG, i_recv_from, MPI_ANY_TAG, MPI_COMM_WORLD, 
					 MPI_STATUS_IGNORE);

			// receive points into new_points
			point_t* new_points = (point_t*) allocate(recv_buff);
			MPI_Recv(new_points, recv_buff, MPI_BYTE, i_recv_from, MPI_ANY_TAG, 
					 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			int num_my_point = sizeof(pt_my_points)/sizeof(point_t);
			int num_new_points = sizeof(new_points)/sizeof(point_t);

			double* d_new_lines = (double*) allocate(sizeof(double)*num_my_point*num_new_points*5);

			int new_line_index = 0;
			for (int j = 0; j < num_my_point; j++){
				for (int k = 0; k < num_new_points; ++k)
				{
					double length = distance(&pt_my_points[k], &new_points[j]);
					d_my_lines[new_line_index + X0] = pt_my_points[j].x;
					d_my_lines[new_line_index + Y0] = pt_my_points[j].y;
					d_my_lines[new_line_index + X1] = new_points[k].x;
					d_my_lines[new_line_index + Y1] = new_points[k].y;
					d_my_lines[new_line_index + LEN] = length;

					new_line_index +=5;
				}
			}
			double* temp = (double*) allocate(sizeof(d_my_lines)+sizeof(d_new_lines));
			int num_my_linr = sizeof(d_my_lines)/sizeof(line_t);
			memcpy(&d_my_lines, &temp, num_my_linr);
			memcpy(&d_new_lines, &temp[num_my_linr], sizeof(d_new_lines)/sizeof(line_t));

			free(d_my_lines);
			d_my_lines = allocate(sizeof(temp));
			memcpy(temp, d_my_lines, sizeof(temp)/sizeof(line_t));
			free(temp);
			free(d_new_lines);
		}
	}// end for
}


void distrib_lines() {
	long num_of_lines;				// number of lines
	int displs[nprocs];				// Used by root only
	double* d_recv_lines = 0; 		// Used by root only
	int* i_recv_counts;				// Used by root only

	if (my_rank == 0) {
		i_recv_counts = allocate (sizeof(int) * nprocs);
	}

	num_of_lines = sizeof(d_my_lines);
	// send the number of lines the receiver should expect
	MPI_Gather(&num_of_lines, 1, MPI_INT, i_recv_counts, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

	if (my_rank==ROOT) {
		displs[0] = 0;
		long total_line_num = i_recv_counts[0];

		// calculate how many total lines are being sent and the displs
        for (int i=1; i < nprocs; i++) {
           total_line_num += i_recv_counts[i];
           displs[i] = displs[i-1] + i_recv_counts[i-1];
        }
		d_recv_lines = (double*) allocate( total_line_num* sizeof(double)*5);        
	}

	MPI_Gatherv(&d_my_lines, num_of_lines, MPI_BYTE, d_recv_lines, i_recv_counts, 
			    displs, MPI_BYTE, ROOT, MPI_COMM_WORLD);

	long l_num_d_lines;
	long l_base ;
	int remainder;
	long *l_send_counts;
	long l_recv_num;
	long *l_displs;

	if(my_rank==ROOT) {
		l_num_d_lines = sizeof(d_recv_lines) / (sizeof(double) * 5);	// Number of lines (5 doubles)
	 	l_base = l_num_d_lines/nprocs;		// Base number of lines to send (lines being 5 doub;es)
	 	remainder = l_num_d_lines%nprocs;	// if there are any remaining lines after the base amount is spilt up
	 	l_send_counts = (long*) allocate(sizeof(long)*nprocs); // Amount of lines (5 doubles) to send to each process 

	 	// Calculate l_send_counts
	 	for (int i = 0; i < nprocs; i++) {
	 		l_send_counts[i] = l_base*5;
	 		if (remainder) {
	 			l_send_counts[i] += 5;
	 			remainder--;	
	 		}
	 	}

		// build displacement array
		l_displs[0] = 0;
		for (int i = 1; i < nprocs; i++) {
			l_displs[i] = l_displs[i-1] + l_send_counts[i-1];
		}	
	}
	//tell processes how many to expect
	MPI_Scatter(l_send_counts, 1, MPI_LONG, &l_recv_num, 1, MPI_LONG, ROOT, MPI_COMM_WORLD);
	
	d_my_lines = (double*) allocate(l_recv_num*sizeof(double)*5);
	//sent lines
	MPI_Scatterv(d_recv_lines, &l_send_counts, &l_displs, MPI_DOUBLE, &d_my_lines, l_recv_num, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
}




////			  //
/// 	Main 	 ///
// 				////

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

	STOP_TIMER(MPIoverhead)

	
	// Root reads in the lines for given file
	if (my_rank==ROOT) {
		read_points(argv);
	}

	// Root scatters the points
	distrib_points();

	
	START_TIMER(generate)

	gen_lines();
	distrib_lines();

	STOP_TIMER(generate)
	
	  //                                      //
	 //  Sort the lines from small to large  //
	//                                      //
	
	START_TIMER(sort)
//	qsort(lines, i_num_lines, sizeof(line_t), compare);
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
