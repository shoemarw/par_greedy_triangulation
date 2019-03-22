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

#define TAG  1		// Used for MPI_sends/revc
#define ROOT 0		// Used for MPI calls and if statements
#define IMAGE

#define X0  0	// first point's x 			| These are used when lines are stored in arrays of 
#define Y0  1	// first point's y 			| doubles rather than structs.
#define X1  2	// second point's x 		| They are used as displacements, for example 
#define Y1  3	// second point's y 		| line_array[i + Y0] would be the i-th line's 
#define LEN 4	// distance the points 		| y-coordinate for first point in the line.

line_t* triang; // Will store the triangulation. Should only be used by proc 0.
long tlines;    // Will keep track of how many lines are in the triangulation at
                // any given time.
bool finished;  // Will be true when the triangulation is complete.
double IMPOSSIBLE_LINE[5] = {0, 0, 0, 0, -1}; // Processes send this during 


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


void double_array_to_struct(double* arr, line_t* new_arr, long size){
	printf("%ld\n", size);
	long index = 0;
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
		print_line(&new_arr[index]);
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
	long bytes_to_expect;
	MPI_Scatter(i_send_count, 1, MPI_INT, &bytes_to_expect, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

	// calculate points a process is responsible for
	my_point_count = bytes_to_expect/sizeof(point_t);

	// Allocate room for the points about to be received
	pt_my_points = (point_t*) allocate((int)(bytes_to_expect));

	// send each process its points
	MPI_Scatterv(points, i_send_count, i_displs_p, MPI_BYTE, pt_my_points, bytes_to_expect,
             MPI_BYTE, ROOT, MPI_COMM_WORLD);
}



void gen_lines() {
	  //				//
	 // Generate lines //
	// 				  //

	// Calculate how many lines that will be created	
	my_line_count = ((my_point_count)*(my_point_count-1))/2;
	// Allocate room for the lines that will be created
	d_my_lines = (double*) allocate(my_line_count * sizeof(double) * 5);


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
			MPI_Send(pt_my_points, my_point_count*sizeof(point_t), MPI_BYTE, i_send_to, TAG, MPI_COMM_WORLD);

//printf("%lf %lf %lf %lf %lf\n", d_my_lines[0], d_my_lines[1], d_my_lines[2], d_my_lines[3], d_my_lines[4]);

			free(pt_my_points);
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
	
					double length = distance(&pt_my_points[j], &pt_new_points[k]);
					d_new_lines[new_line_index + X0] = pt_my_points[j].x;
					d_new_lines[new_line_index + Y0] = pt_my_points[j].y;
					d_new_lines[new_line_index + X1] = pt_new_points[k].x;
					d_new_lines[new_line_index + Y1] = pt_new_points[k].y;
					d_new_lines[new_line_index + LEN] = length;
// printf("Proc0 'new' lines: ");
// printf("%lf  ", d_new_lines[new_line_index + X0]);
// printf("%lf  ", d_new_lines[new_line_index + Y0]);
// printf("%lf  ", d_new_lines[new_line_index + X1]);
// printf("%lf  ", d_new_lines[new_line_index + Y1]);
// printf("%lf \n", d_new_lines[new_line_index + LEN]);

					new_line_index +=5;
				}
			} // end out for


			// merge d_my_lines with d_new_lines
			double* temp_doubles = (double*) allocate((my_line_count+new_line_count)*sizeof(double)*5);
			for (int i = 0; i < (my_line_count*5); i++) {
				temp_doubles[i] = d_my_lines[i];
			}
			for(int i = 0; i < (new_line_count*5); i++) {
				temp_doubles[i + my_line_count*5] = d_new_lines[i];
			}
			free(d_my_lines);
			free(d_new_lines);
			d_my_lines = temp_doubles;

// //print the stuff in d_my_lines
// for (int i = 0; i < (my_line_count+new_line_count); i++) { 
// //printf("Line p00p %lf  %lf  %lf  %lf  %lf  \n", d_my_lines[0+5*i],d_my_lines[1+5*i],
// d_my_lines[2+5*i],d_my_lines[3+5*i],d_my_lines[4+5*i]);
// }

			// merge pt_my_points with pt_new_points
			point_t* temp_points = (point_t*) allocate((my_point_count+point_recv_count)*sizeof(point_t));
			for (int i = 0; i < my_point_count; i++) {
				temp_points[i] = pt_my_points[i];
			}
			for(int i = 0; i < point_recv_count; i++) {
				temp_points[i + my_point_count] = pt_new_points[i];
			}
			free(pt_my_points);
			free(pt_new_points);
			pt_my_points = temp_points;

			// update my_line_count	
			my_line_count += new_line_count;

for(int i = 0; i < my_line_count; i++) {

}

		} // end of receiver branch of if
	}// end for
}// end of gen_lines


void distrib_lines() {
	int displs[nprocs];		// Used by root only
	double* d_recv_lines;	// Used by root only
	int* i_recv_counts;		// Used by root only

	if (my_rank == ROOT) {
		i_recv_counts = allocate (sizeof(int) * nprocs);
	}

	// send the number of lines a process will be sending on the gatherv
	MPI_Gather(&my_line_count, 1, MPI_INT, i_recv_counts, 1, MPI_INT, ROOT, MPI_COMM_WORLD);
	long total_line_num;
	if (my_rank == ROOT) {
		for(int i = 0; i < nprocs; i++) {
			i_recv_counts[i] *= 5;
		}

		displs[0] = 0;
		total_line_num = i_recv_counts[0];

		// calculate how many total lines are being sent and the displs
        for (int i=1; i < nprocs; i++) {
           	total_line_num += i_recv_counts[i];
           	displs[i] = displs[i-1] + i_recv_counts[i-1];
        }
		d_recv_lines = (double*) allocate(total_line_num* sizeof(double));

// printf("Proc0 is expecting %d doubles from proc0\n", i_recv_counts[0]);
// printf("Proc0 is expecting %d doubles from proc1\n", i_recv_counts[1]);
// printf("Proc0 is making room for %ld doubles\n", total_line_num);

	}
// printf("I am prco %d and I am sending %ld doubles\n", my_rank, my_line_count*5);
// printf("I am prco %d and I am sending:\n", my_rank);
// printf("Line %lf  %lf  %lf  %lf  %lf  \n", d_my_lines[0+5],d_my_lines[1+5],d_my_lines[2+5],d_my_lines[3+5],d_my_lines[4+5]);

	MPI_Gatherv(d_my_lines, (my_line_count*5), MPI_DOUBLE, d_recv_lines, i_recv_counts, 
			    displs, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

// if (my_rank==ROOT){
// 	for (int i = 0; i < total_line_num; i+=5)
// 	{
// 		printf("Line: %d\n", (i/5)+1);
// 		for (int j = 0; j < 5; ++j)
// 		{
// 			printf("%lf  ", d_recv_lines[i+j]);
// 		}
// 		printf("\n");
// 	}
// }



	long l_base;
	int remainder;
	int *i_send_counts;
	int i_recv_doubs;
	int *i_displs;

	if(my_rank==ROOT) {
		free(i_recv_counts);
		my_line_count = total_line_num/5;

	 	l_base = my_line_count/nprocs;		// Base number of lines to send (lines being 5 doubles)
	 	remainder = my_line_count%nprocs;	// if there are any remaining lines after the base amount is split up
	 	i_send_counts = (int*) allocate(sizeof(int) * nprocs); // Amount of lines (5 doubles) to send to each process 
	 	// Calculate i_send_counts
	 	for (int i = 0; i < nprocs; i++) {
	 		i_send_counts[i] = l_base*5;	// (*5) is to account for lines being five doubles
	 		if (remainder) {
	 			i_send_counts[i] += 5; 	// +5 because each line is really 5 doubles at this point
	 			remainder--;	
	 		}
	 	}


		// build displacement array
		i_displs = (int *) allocate(sizeof(int)*nprocs);
		i_displs[0] = 0;
		for (int i = 1; i < nprocs; i++) {
			i_displs[i] = i_displs[i-1] + i_send_counts[i-1];
		}
	}
	// tell processes how many doubles to expect in the scatterv
	MPI_Scatter(i_send_counts, 1, MPI_INT, &i_recv_doubs, 1, MPI_INT, ROOT, MPI_COMM_WORLD);

	// calculate how many lines the process is responsible for
	my_line_count = i_recv_doubs/5;

	// create room for for the lines
	d_my_lines = (double*) allocate(i_recv_doubs*sizeof(double));

	// scatter lines
	MPI_Scatterv(d_recv_lines, i_send_counts, i_displs, MPI_DOUBLE, d_my_lines, i_recv_doubs, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

	ln_my_lines = (line_t *) allocate(my_line_count*sizeof(line_t));

	double_array_to_struct(d_my_lines, ln_my_lines, my_line_count);

// if(my_rank==ROOT) {
// 	printf("after double_array_to_struct\n");
// 	for (int i = 0; i < my_line_count; i++) {
// 		printf("I am proc %d and this is line %d\n", my_rank, i);
// 		print_line(&ln_my_lines[i]);
// 	}
// }
}


/*
 * Builds the Greedy Triangulation. 
 */
void triangulate() {
	// The triangulation will be stored as an array of lines. The triangulation
	// is built on process zero iteratively as successive global minimal lines
	// are found. Allocate enough space to potentially hold every line.
	triang = (line_t*) allocate(my_line_count*sizeof(line_t));

	// my_unknown is each local processes' number of lines whose status is 
	// unknown. So it counts the number of lines that may or may not belong
	// to the triangulation that a give process has. Initially all of a 
	// processes' lines have unknown status.
	long my_unknown = my_line_count;
	tlines = 0; // Set the number of lines currently in the triangulation to 0.
	// Keep participating in global communications until all processes have
	// resolved the status of their local set of lines.
	while (!finished) {
		// If this process still has lines of unknown status it must
		// work to resolve them.
		if (my_unknown > 0) {
			// Convert this processes' minimal (smallest) line to an array of
			// five doubles for Allgather.
			double my_min_line[5];
			break;
			point_t p = *(ln_my_lines[0].p);   /// Make sure there is a line in ln_my_lines
			point_t q = *(ln_my_lines[0].q);
			
			my_min_line[0] = p.x;
			my_min_line[1] = p.y;
			my_min_line[2] = q.x;
			my_min_line[3] = q.y;
			// Prepare an array to receive each processe's minimal line.
			double* recv_buf = (double*) allocate(5*nprocs);                    ///// Should be 5*nprocs*sizeof(double) ////
			// Make sure each process has an array of each processes' min line.
			MPI_Allgather(my_min_line, 5, MPI_DOUBLE, 
				          recv_buf, 5, MPI_DOUBLE, MPI_COMM_WORLD);
			// Find the global minimal line.
			int min_line_index = 0; // Will hold index of the global min line.
			for (int i = 0; i < nprocs; i++) {
				// Compare the length of the current smallest line to the
				// i^th line's length. If the length is not positive ignore it
				// because it was a special value sent from a process with no
				// more lines of unknown status.
				if ((recv_buf[i*4] > 0) && 
					(recv_buf[i*4] < recv_buf[min_line_index*4])) {
					min_line_index = i;
				}
			}
			// Will hold the minimal line.
			line_t* min_line;
			// Lets us know whether or not this processes min line was
			// included in the triangulation or not.
			int start;
			// See if this process has the global min line, if so we must
			// adjust its number of lines of unknown status and set min_line.
			if (my_rank == min_line_index) {
				my_unknown--;
				min_line = &ln_my_lines[0];
				start = 1; // This processes' min was used.
			// Otherwise we must build the min_line from data in the recv_buf
			// (While making sure to use the appropriate index!)
			} else {
				// Get the minimal line
				min_line = (line_t*) allocate(sizeof(line_t));
				point_t *p = (point_t*) allocate(sizeof(point_t));
				point_t *q = (point_t*) allocate(sizeof(point_t));
				p->x = recv_buf[min_line_index*4];
				p->y = recv_buf[min_line_index*4 + 1];
				q->x = recv_buf[min_line_index*4 + 2];
				q->y = recv_buf[min_line_index*4 + 3];
				min_line->p = p;
				min_line->q = q;
				free(p);
				free(q);
				start = 0; // This processes' min was not used.
			}
			// Have process zero add min_line to the triangulation.
			if (my_rank == 0) {
				triang[tlines] = *min_line;
			}
			// Free the receive buffer
			free(recv_buf);
			// Allocate an array of lines to hold the lines that don't 
			// intersect with the global minimum.
			line_t* temp = (line_t*) allocate(my_line_count*sizeof(line_t));
			int end = my_unknown;
			int temp_size = 0;

			for (int j = start; j < end; j++) {
				// Run the intersection test and only include lines which dont
				// conflict with the global min (min_line) It is ok if a line
				// shares endpoints with min_line
				if (share_endpoint(min_line, &ln_my_lines[j]) ||
					 !intersects(min_line, &ln_my_lines[j])) {
					temp[temp_size] = ln_my_lines[j];
					temp_size++;
				} else {
					my_unknown--;
				}	
			}
			// Write all of the valid lines from temp to ln_my_lines to prepare
			// for the next iteration.		
			copy_array(temp, ln_my_lines, temp_size);
			free(temp);
		// If this process has no more lines of unknown status then it must still
	    // participate in global communications to avoid deadlock. Have it send
		// a special value to the other processes (which they will ignore). The 
		// special value is a line whose endpoints are at the origin and it has a
		// distance of -1.
		} // end if (my_unknown > 0)
		else {
			// Prepare an array to receive each processe's minimal line.
			double* recv_buf = (double*) allocate(5*nprocs);
			MPI_Allgather(IMPOSSIBLE_LINE, 5, MPI_DOUBLE, 
				          recv_buf, 5, MPI_DOUBLE, MPI_COMM_WORLD);
			// Check if all of the lines have distance of -1. If so then we are
			// done!
			int count = 0;
			for (int i = 0; i < nprocs; i++) {
				if (recv_buf[i*4] == -1) {
					count++;
				}
			}
			if (count == nprocs) {
				finished = true;
			}
			free(recv_buf); 
		} // end else
	} // end while
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
	
	// Root reads in the lines for given file
	if (my_rank==ROOT) {
		read_points(argv);
	}

	// Root scatters the points
	distrib_points();

	MPI_Barrier(MPI_COMM_WORLD);
	START_TIMER(generate)
	gen_lines();
	distrib_lines();
	MPI_Barrier(MPI_COMM_WORLD);
	STOP_TIMER(generate)

// if(my_rank==ROOT) {
// 	printf("Before sort\n");
// 	for (int i = 0; i < my_line_count; i++) {
// 		printf("I am proc %d and this is line %d\n", my_rank, i);
// 		print_line(&ln_my_lines[i]);
// 	}
// }
	  //                                      //
	 //  Sort the lines from small to large  //
	//                                      //
	MPI_Barrier(MPI_COMM_WORLD);
	START_TIMER(sort)
	qsort(ln_my_lines, (sizeof(ln_my_lines)/sizeof(line_t)), sizeof(line_t), compare);
	MPI_Barrier(MPI_COMM_WORLD);
	STOP_TIMER(sort)
// if(my_rank==ROOT) {
// 	for (int i = 0; i < my_line_count; i++) {
// 		printf("I am proc %d and this is line %d\n", my_rank, i);
// 		print_line(&ln_my_lines[i]);
// 	}
// }
	
	  //                                   //
     //  Greedily build the tringulation  //
    //	                                 //
	
	MPI_Barrier(MPI_COMM_WORLD);
	START_TIMER(triangulate)
	triangulate();

	MPI_Barrier(MPI_COMM_WORLD);
	STOP_TIMER(triangulate)
	
	   //                                        //
      //  Triangulation Done, Display Stats     //
     //	  Create image & Output File in proc0  //
	//                                        //

	// These stats are only for the portions of the code specific to the three
	// phases of building the greedy triangulation. Generate all lines, sort the
	// lines in non-decreasing order, and greedily adding line segments to the
	// triangulation.
	if (my_rank==ROOT) {
		printf("Gent: %.4f  Sort: %.4f  Tria: %.4f \n",
	        GET_TIMER(generate), GET_TIMER(sort), 
	        GET_TIMER(triangulate));

# ifdef IMAGE	
		generate_image(triang, tlines);
# endif

		// Store the triangulation in a file. Each line in the file corresponds to 
		// a line in the triangulation. Each line consists of two tuples representing
		// the end points of the line.
		FILE* write_file = open_file("triangle_result.txt", "w");
		
		// The first line of the file specifies the number of lines in the file.
		fprintf(write_file, "%ld\n", tlines);
		
		// Write triangulation to file for turtle processing
		for (int i = 0; i < tlines; i++)
		{
		    point_t p = *(triang[i].p);
		    point_t q = *(triang[i].q);
		    fprintf(write_file, "(%lf, %lf) (%lf, %lf)\n", p.x, p.y, q.x, q.y);
		}
		
		fclose(write_file);
	}
	
	// Clean up and exit
	if (my_rank == ROOT)
		free(points);

	free(triang);
	MPI_Finalize();
	return (EXIT_SUCCESS);
}
