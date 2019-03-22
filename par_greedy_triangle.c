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

const int  TAG = 1;		// Used for MPI_sends/revc
const int  ROOT = 0;	// Used for MPI calls and if statements

const int X0 = 0;	// first point's x 			| These are used when lines are stored in arrays of 
const int Y0 = 1;	// first point's y 			| doubles rather than structs.
const int X1 = 2;	// second point's x 		| They are used as displacements, for example 
const int Y1 = 3;	// second point's y 		| line_array[i + Y0] would be the i-th line's 
const int LEN = 4;	// distance the points 		| y-coordinate for first point in the line.

int nprocs;         	// number of processes
int my_rank;        	// rank of a particular process
point_t* pt_my_points; 	// a processes' portion of the point set.
line_t*  ln_my_lines;  	// a processes' portion of lines.
double* d_my_lines;		// a processes' portion of lines while in double format.


long my_point_count;	// Count of how many points a processes is responsible for
long my_line_count = 0; // Count of how many lines a processes is responsible for

//TEMPORARILY MAKING THESE GLOBAL, THIS MAY NEED TO BE CHANGED
long l_num_points;
point_t* points;
line_t* lines;


void double_array_to_struct(double* arr, line_t* new_arr, long size){
	long index = 0;
if (my_rank==ROOT) printf("%ld\n", size);
	for (long i = 0; i < size; i+=5) {
		point_t *p0 = (point_t*) allocate(sizeof(point_t));
		point_t *p1 = (point_t*) allocate(sizeof(point_t));
		p0->x = arr[i+X0];
		p0->y = arr[i+Y0];
		p1->x = arr[i+X1];
		p1->y = arr[i+Y1];

		line_t* l = (line_t*) allocate(sizeof(line_t));
		// set the values of the line and store it.
		l->p = p0;
		l->q = p1;
		l->len = arr[i+LEN];
		new_arr[index] = *l;
if (my_rank==ROOT) printf("hello from line 65\n");
if (my_rank==ROOT) print_line(l);
		index++;
		free(l);
	}
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
	int i_send_count[nprocs];	// an array of how bytes each process will receive; significant only to root
	int i_displs_p[nprocs];			// the displacements for the scatterv; significant only to root

	if (my_rank==ROOT) {
		// use integer division to determine the base amount for points each process will receive 
		long base_point_count = l_num_points/(long)nprocs;

		// get the remainder to see how many leftover points there are
		int remainder = l_num_points%nprocs;

		// fill the array with the base number, then if there are remainders left add one to the 
		// count of how many points the process will receive.
		for (int i = 0; i < nprocs; i++) {
			i_send_count[i] = base_point_count * sizeof(point_t);
			if (remainder > 0) {
				i_send_count[i] += 1 * sizeof(point_t);
				remainder--;
			} // end if
		} // end for

		// build displacement array
		i_displs_p[0] = 0;
		for (int i = 1; i < nprocs; i++) {
			i_displs_p[i] = i_displs_p[i-1] + i_send_count[i-1];
		}
	}
	// send each process how many points it should expect
	MPI_Scatter(i_send_count, 1, MPI_INT, &my_point_count, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

	// Allocate room for the points about to be recieved
	pt_my_points = (point_t*) allocate((int)(my_point_count*sizeof(point_t)));

	// send each process its points
	MPI_Scatterv(points, i_send_count, i_displs_p, MPI_BYTE, pt_my_points, my_point_count,
             MPI_BYTE, ROOT, MPI_COMM_WORLD);
}



void gen_lines() {
	  //				//
	 // Generate lines //
	// 				  //

	// Calculate how many lines that will be created	
	int i_num_lines = ((my_point_count-1)*(my_point_count-2))/2;

	// Allocate room for the lines that will be created
	d_my_lines = (double*) allocate(i_num_lines * sizeof(double) * 5);


	long l_index = 0;

	// Create lines in double format, i.e. a line is 5 doubles.
	for (int i = 0; i < my_point_count; i++) {
		// Compute the distance between point i and every point
		// from i+1 onward. Then 'make' and store the corresponding
		// line.
		for (int j = i+1; j < my_point_count; j++) {
			double length = distance(&pt_my_points[i], &pt_my_points[j]);
			d_my_lines[l_index + X0] = pt_my_points[i].x;
			d_my_lines[l_index + Y0] = pt_my_points[i].y;
			d_my_lines[l_index + X1] = pt_my_points[j].x;
			d_my_lines[l_index + Y1] = pt_my_points[j].y;
			d_my_lines[l_index + LEN] = length;	

			// Increment index to next line
			l_index +=5;
		}
	}
	  //													  //
	 // Start sending/receiving points to create more lines  //
	//														//

	int point_recv_count; 	// Used to check how many objects will be sent in MPI_send
	double* d_new_lines;
	point_t* pt_new_points;
	
	// In this for loop we calculated all the remaining lines - which means we'll
	// need to send/receive points from other processes. This is implemented
	// in a binomial tree structure.
	for (int iteration_square = 1; iteration_square < nprocs; iteration_square *= 2) {
		// If process is a sender this iteration:
		if (my_rank&iteration_square) {

			// calculate the process number of whom to send to
			int i_send_to = (my_rank-iteration_square+nprocs)%nprocs;

			// send the number of points the receiver should expect
			MPI_Send(&my_point_count, 1, MPI_LONG, i_send_to, TAG, MPI_COMM_WORLD);

			// send the points
			MPI_Send(pt_my_points, my_point_count, MPI_BYTE, i_send_to, TAG, MPI_COMM_WORLD);
			
			break;  // done, nothing left for this process to do in this function
		}
		// If process is a receiver this iteration:
		else {

			// calculate the process number of whom to receive from
			int i_recv_from = my_rank+iteration_square; // The process number to receive from
			
			// receive the number of points about to get sent
			MPI_Recv(&point_recv_count, 1, MPI_LONG, i_recv_from, MPI_ANY_TAG, MPI_COMM_WORLD, 
					 MPI_STATUS_IGNORE);

			// calculate how many points that will be received 
			long bytes_to_recv = point_recv_count*sizeof(point_t);

			// create a array for bytes/points about to be received
			pt_new_points = (point_t*) allocate(bytes_to_recv);

			// receive points into pt_new_points
			MPI_Recv(pt_new_points, bytes_to_recv, MPI_BYTE, i_recv_from, MPI_ANY_TAG, 
					 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			// calculate number of new lines to be created
			long new_line_count = my_point_count*point_recv_count;

			// create array to store lines about to be created
			d_new_lines = (double*) allocate(sizeof(double)*new_line_count*5);

			// create long to keep tract of d_new_lines index through for loop
			long new_line_index = 0;

			// create new lines by taking one of the process's current lines and making lines from that point
			// to all of the newly received points. This is done in a double for loop.
			for (int j = 0; j < my_point_count; j++){
				for (int k = 0; k < point_recv_count; ++k)
				{
					double length = distance(&pt_my_points[k], &pt_new_points[j]);
					d_my_lines[new_line_index + X0] = pt_my_points[j].x;
					d_my_lines[new_line_index + Y0] = pt_my_points[j].y;
					d_my_lines[new_line_index + X1] = pt_new_points[k].x;
					d_my_lines[new_line_index + Y1] = pt_new_points[k].y;
					d_my_lines[new_line_index + LEN] = length;

					new_line_index +=5;
				}
			} // end out for

			// bunch of ugly code to combine both arrays of points into pt_my_points /////
			point_t *pt_temp = (point_t *) allocate(my_point_count + point_recv_count);	//
			memcpy(pt_my_points, pt_temp, my_point_count);								//
			memcpy(pt_new_points, &pt_temp[my_point_count], point_recv_count);			//
			pt_my_points = (point_t *) allocate(my_point_count + point_recv_count);		//
			memcpy(pt_temp,pt_my_points, (my_point_count + point_recv_count));			//
			free(pt_temp);																//
			// free(pt_new_points);														//
			//////////////////////////////////////////////////////////////////////////////

			// update my_point_count
			my_point_count += point_recv_count;

			// bunch of ugly code to combine both arrays of lines into d_my_lines ////////
			double* d_temp = (double*) allocate(my_line_count+new_line_count);			//
			memcpy(&d_my_lines, &d_temp, my_line_count);								//
			memcpy(&d_new_lines, &d_temp[my_line_count], my_line_count);				//
			d_my_lines = allocate(my_line_count+new_line_count);						//
			memcpy(d_temp, d_my_lines, (my_line_count+new_line_count));					//
			free(d_temp);																//
			// free(d_new_lines);														//
			//////////////////////////////////////////////////////////////////////////////

			// update my_line_count
			my_line_count += new_line_count;

		} // end of receiver branch of if
	}// end for
}// end of gen_lines


void distrib_lines() {
	int displs[nprocs];				// Used by root only
	double* d_recv_lines = 0; 		// Used by root only
	int* i_recv_counts;				// Used by root only

	if (my_rank == 0) {
		i_recv_counts = allocate (sizeof(int) * nprocs);
	}

	// send the number of lines a process will be sending on the gatherv
	MPI_Gather(&my_line_count, 1, MPI_INT, i_recv_counts, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

	if (my_rank==ROOT) {
		displs[0] = 0;
		long total_line_num = i_recv_counts[0];

		// calculate how many total lines are being sent and the displs
        for (int i=1; i < nprocs; i++) {
           total_line_num += i_recv_counts[i];
           displs[i] = displs[i-1] + i_recv_counts[i-1];
        }
		d_recv_lines = (double*) allocate(total_line_num* sizeof(double)*5);        
	}

	MPI_Gatherv(&d_my_lines, my_line_count, MPI_BYTE, d_recv_lines, i_recv_counts, 
			    displs, MPI_BYTE, ROOT, MPI_COMM_WORLD);

	long l_base;
	int remainder;
	long *l_send_count;
	long l_recv_doubs;
	int *i_displs;

	if(my_rank==ROOT) {
		my_line_count = 0;
		for (int i = 0; i < nprocs; i++) {
			my_line_count += d_recv_lines[i];
		}


	 	l_base = my_line_count/nprocs;		// Base number of lines to send (lines being 5 doubles)
	 	remainder = my_line_count%nprocs;	// if there are any remaining lines after the base amount is split up
	 	l_send_count = (long*) allocate(sizeof(long) * nprocs); // Amount of lines (5 doubles) to send to each process 

	 	// Calculate l_send_count
	 	for (int i = 0; i < nprocs; i++) {
	 		l_send_count[i] = l_base*5;	// (*5) is to account for lines being five doubles
	 		if (remainder) {
	 			l_send_count[i] += 5; 	// +5 because each line is really 5 doubles at this point
	 			remainder--;	
	 		}
	 	}

		// build displacement array
		i_displs = (int *) allocate(sizeof(int)*nprocs);
		i_displs[0] = 0;
		for (int i = 1; i < nprocs; i++) {
			i_displs[i] = i_displs[i-1] + l_send_count[i-1];
		}	
	}
	// tell processes how many lines to expect in the scatterv
	MPI_Scatter(l_send_count, 1, MPI_LONG, &l_recv_doubs, 1, MPI_LONG, ROOT, MPI_COMM_WORLD);
	
	// calculate how many lines the process is responsible for
	my_line_count = l_recv_doubs/5;

	// create room for for the lines
	d_my_lines = (double*) allocate(l_recv_doubs*sizeof(double));

	// scatter lines
	MPI_Scatterv(d_recv_lines, (int *)l_send_count, i_displs, MPI_DOUBLE, &d_my_lines, l_recv_doubs, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

	ln_my_lines = (line_t *) allocate(my_line_count*sizeof(line_t));

	double_array_to_struct(d_my_lines, ln_my_lines, my_line_count);
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


	// if (my_rank==0){	
	// 	for (int i = 0; i < 10; i++) {
	// 		print_line(&ln_my_lines[i]);
	// 	}
	// }




	
	  //                                      //
	 //  Sort the lines from small to large  //
	//                                      //
	
	START_TIMER(sort)
	// qsort(ln_my_lines, (sizeof(ln_my_lines)/sizeof(line_t)), sizeof(line_t), compare);
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
		printf("Gent: %.4f  Sort: %.4f  Tria: %.4f Overhead: %.4f\n",
	        GET_TIMER(generate), GET_TIMER(sort), 
	        GET_TIMER(triangulate), GET_TIMER(MPIoverhead));
	}
	
	// Clean up and exit
	// free(points);
	// free(lines);
	MPI_Finalize();
	return (EXIT_SUCCESS);
}
