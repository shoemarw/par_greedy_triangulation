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

////////////////////////////////////////////////////////////////////////////////
line_t* triang; // Will store the triangulation. Should only be used by proc 0.
long tlines;    // Will keep track of how many lines are in the triangulation at
                // any given time.
bool finished;  // Will be true when the triangulation is complete.
double[5] IMPOSSIBLE_LINE = [0, 0, 0, 0, -1]; // Processes send this during 
                                           // Allgather if they are finished.
////////////////////////////////////////////////////////////////////////////////


void double_array_to_struct(double* arr, line_t* new_arr, long size){
	const int X0 = 0;
	const int Y0 = 1;
	const int X1 = 2;
	const int Y1 = 3;
	const int LEN = 4;
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

	int i_send_counts[nprocs];	// an array of how many points each process will recieve / how many root sends
	int i_displs_p[nprocs];		// the displacements for the scatterv, significant only to root
	int i_points_to_recv = 0;	// used to store a single number from i_send_counts

	if (my_rank==ROOT) {
		// use interger division to determine the base amount for points each process will recieve 
		long base_point_count = l_num_points/(long)nprocs;
		printf("l_num_points %ld\n", l_num_points);

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
			printf("i_send_counts %d\n", i_send_counts[i]);
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
	int i_num_lines = ((l_num_points-1)*(l_num_points-2))/2;
	d_my_lines = (double*) allocate(i_num_lines * sizeof(double) * 5);
if(my_rank==ROOT) {printf("166 %ld\n", sizeof(*d_my_lines));}
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

			// receive points into pt_new_points
			point_t* pt_new_points = (point_t*) allocate(recv_buff);
			MPI_Recv(pt_new_points, recv_buff, MPI_BYTE, i_recv_from, MPI_ANY_TAG, 
					 MPI_COMM_WORLD, MPI_STATUS_IGNORE);

			int num_my_point = sizeof(pt_my_points)/sizeof(point_t);
			int num_pt_new_points = sizeof(pt_new_points)/sizeof(point_t);

			double* d_new_lines = (double*) allocate(sizeof(double)*num_my_point*num_pt_new_points*5);

			int new_line_index = 0;
			for (int j = 0; j < num_my_point; j++){
				for (int k = 0; k < num_pt_new_points; ++k)
				{
					double length = distance(&pt_my_points[k], &pt_new_points[j]);
					d_my_lines[new_line_index + X0] = pt_my_points[j].x;
					d_my_lines[new_line_index + Y0] = pt_my_points[j].y;
					d_my_lines[new_line_index + X1] = pt_new_points[k].x;
					d_my_lines[new_line_index + Y1] = pt_new_points[k].y;
					d_my_lines[new_line_index + LEN] = length;

					new_line_index +=5;
				}
			}
if(my_rank==ROOT) {printf("235 %ld\n", sizeof(d_my_lines));}
			long l_size_my_points = sizeof(pt_my_points)/sizeof(point_t);
			long l_size_new_points = sizeof(pt_new_points)/sizeof(point_t);
			point_t *pt_temp = (point_t *) allocate(sizeof(pt_my_points) + sizeof(pt_new_points));
			memcpy(pt_my_points, pt_temp, l_size_my_points);
			memcpy(pt_new_points, &pt_temp[l_size_my_points],l_size_new_points);
			free(pt_new_points);
			free(pt_my_points);
			pt_my_points = (point_t *) allocate(sizeof(pt_temp));
			memcpy(pt_temp,pt_my_points,sizeof(pt_my_points)/sizeof(point_t));
			free(pt_temp);

			double* d_temp = (double*) allocate(sizeof(d_my_lines)+sizeof(d_new_lines));
			int num_my_linr = sizeof(d_my_lines)/sizeof(line_t);
			memcpy(&d_my_lines, &d_temp, num_my_linr);
			memcpy(&d_new_lines, &d_temp[num_my_linr], sizeof(d_new_lines)/sizeof(line_t));

			free(d_my_lines);
			d_my_lines = allocate(sizeof(d_temp));
			memcpy(d_temp, d_my_lines, sizeof(d_temp)/sizeof(line_t));
			free(d_temp);
			free(d_new_lines);
		}
	}// end for
if(my_rank==ROOT) {printf("259 %ld\n", sizeof(d_my_lines));}
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
	int *i_send_counts;
	long l_recv_num;
	int *i_displs;

	if(my_rank==ROOT) {
		l_num_d_lines = sizeof(d_recv_lines) / (sizeof(double) * 5);	// Number of lines (5 doubles)
	 	l_base = l_num_d_lines/nprocs;		// Base number of lines to send (lines being 5 doub;es)
	 	remainder = l_num_d_lines%nprocs;	// if there are any remaining lines after the base amount is spilt up
	 	i_send_counts = (int*) allocate(sizeof(int) * nprocs); // Amount of lines (5 doubles) to send to each process 

	 	// Calculate i_send_counts
	 	for (int i = 0; i < nprocs; i++) {
	 		i_send_counts[i] = l_base*5;
	 		if (remainder) {
	 			i_send_counts[i] += 5;
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
	//tell processes how many to expect
	MPI_Scatter(i_send_counts, 1, MPI_LONG, &l_recv_num, 1, MPI_LONG, ROOT, MPI_COMM_WORLD);
	
	d_my_lines = (double*) allocate(l_recv_num*sizeof(double)*5);
	//sent lines
	MPI_Scatterv(d_recv_lines, i_send_counts, i_displs, MPI_DOUBLE, &d_my_lines, l_recv_num, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

	ln_my_lines = (line_t *) allocate((l_recv_num/5)*sizeof(line_t));
if (my_rank==ROOT) {
printf("hello from line 326\n");
printf("%ld\n", sizeof(d_my_lines));
printf("%ld\n", sizeof(double));
}

	double_array_to_struct(d_my_lines, ln_my_lines, sizeof(d_my_lines)/(sizeof(double) * 5));
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
	qsort(ln_my_lines, (sizeof(ln_my_lines)/sizeof(line_t)), 
		                sizeof(line_t), compare);
	STOP_TIMER(sort)
	
	  //                                   //
     //  Greedily build the tringulation  //
    //	                                 //
	
	START_TIMER(triangulate)

	// The triangulation will be stored as an array of lines. The triangulation
	// is built on process zero iteratively as successive global minimal lines
	// are found. Allocate enough space to potentially hold every line.
	triang = (line_t*) allocate(num_line_count*sizeof(line_t));

	// my_unknown is each local processes' number of lines whose status is 
	// unknown. So it counts the number of lines that may or may not belong
	// to the triangulation that a give process has. Initially all of a 
	// processes' lines have unknown status.
	long my_unknown = my_line_count;
	tlines = 0; // Set the number of lines currently in the triangulation to 0.
	// Keep participating in global communications until all processes have
	// resolved the status of their local set of lines.
	while (!finished) {
		// If this process still has lines of unkown status it must
		// work to resolve them.
		if (my_unknown > 0) {
			// Convert this processes' minimal (smallest) line to an array of
			// five doubles for Allgather.
			double my_min_line[5];
			point_t p = *(ln_my_lines[0].p);
			point_t q = *(ln_my_lines[0].q);
			my_min_line[0] = p.x;
			my_min_line[1] = p.y;
			my_min_line[2] = q.x;
			my_min_line[3] = q.y;
			// Prepare an array to receive each processe's minimal line.
			double* recv_buf = (double*) allocate(5*nprocs);
			// Make sure each process has an array of each processes' min line.
			MPI_Allgather(my_min_line, 5, MPI_DOUBLE, 
				          recv_buf, 5, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
			// Find the global minimal line.
			int min_line_index = 0; // Will hold index of the global min line.
			for (int i = 0; i < nprocs, i++) {
				// Compare the lenth of the current smallest line to the
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
				min_line = ln_my_lines[0];
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
				triang[tlines] = min_line;
			}
			// Free the receive buffer
			free(recv_buf);
			// Allocate an array of lines to hold the lines that don't 
			// intersect with the global minimum.
			line_t* temp = (line_t*) allocate(my_line_count*sizeof(line_t));
			int end = my_unknown;
			int temp_size = 0;

			for (int j = start; j < end, j++) {
				// Run the intersection test and only include lines which dont
				// conflict with the global min (min_line) It is ok if a line
				// shares endpoints with min_line
				if (share_endpoint(min_line, ln_my_lines[j]) ||
					 !intersects(min_line, ln_my_lines[j])) {
					temp[temp_size] = ln_my_lines[j]
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
		} else {
			// Prepare an array to receive each processe's minimal line.
			double* recv_buf = (double*) allocate(5*nprocs);
			MPI_Allgather(IMPOSSIBLE_LINE, 5, MPI_DOUBLE, 
				          recv_buf, 5, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);
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
		}
	}
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

	// In process zero produce the image if IMAGE is defined

	// In process zero write the triangulation to an output file
	
	// Clean up and exit
	// free(points);
	// free(lines);
	free(triang);
	MPI_Finalize();
	return (EXIT_SUCCESS);
}