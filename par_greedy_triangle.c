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
bool finished = false;  // Will be true when the triangulation is complete.
double IMPOSSIBLE_LINE[5] = {0, 0, 0, 0, -1}; // Processes send this during 


int nprocs;         	// number of processes
int my_rank;        	// rank of a particular process
line_t* ln_my_lines;  	// a processes' portion of lines.
double* d_my_lines;		// a processes' portion of lines while in double format.
double* d_all_lines;	// all lines in double format, used by root

long my_line_count = 0; // Count of how many lines a processes is responsible for

//TEMPORARILY MAKING THESE GLOBAL, THIS MAY NEED TO BE CHANGED
long l_num_points;
point_t* points;
point_t* points_to_free;
point_t* min_line_points;
long mlp_index = 0;			// index of next point to be stored


void double_array_to_struct(double* arr, line_t* new_arr, long size){
	points_to_free = (point_t*) allocate(sizeof(point_t)*2*size);
	int points_index = 0;

	long index = 0;
	for (long i = 0; i < size*5; i+=5) {
		line_t* l = (line_t*) allocate(sizeof(line_t));
		
		// Convert first point
		point_t *p0 = (point_t*) allocate(sizeof(point_t));
		p0->x = arr[i+X0];
		p0->y = arr[i+Y0];

		// store point in array to free later
		points_to_free[points_index] = *p0;

		// Put first point into struct
		l->p = &points_to_free[points_index];
		free(p0);
		points_index++;


		// Convert second point
		point_t *p1 = (point_t*) allocate(sizeof(point_t));
		p1->x = arr[i+X1];
		p1->y = arr[i+Y1];

		// store point in array to free later
		points_to_free[points_index] = *p1;

		// Put second point into struct
		l->q = &points_to_free[points_index];
		free(p1);
		points_index++;	


		// Put length into struct
		l->len = arr[i+LEN];
		
		new_arr[index] = *l;
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
	long i = 0;      // Index for storing points.

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

void gen_lines() {
	  //                      //
	 //  Generate all lines  //
	//                      //
	
	// Make all possible line segments between the points
	// and compute the length of each line.
	// Compute the number of lines, note that the below formula
	// will always resolve as an int because one of l_num_points or
	// l_num_points-1 will be divisible by 2.
	my_line_count = ((l_num_points)*(l_num_points-1))/2;
	d_all_lines = (double*) allocate(sizeof(double)*my_line_count*5);	

	long index = 0;
	for (int i = 0; i < l_num_points; i++) {
		// Compute the distance between point i and every point
		// from i+1 onward. Then 'make' and store the corresponding
		// line.
		for (int j = i+1; j < l_num_points; j++) {

			double length = distance(&points[i], &points[j]);

			d_all_lines[index*5+X0] = points[i].x;
			d_all_lines[index*5+Y0] = points[i].y;
			d_all_lines[index*5+X1] = points[j].x;
			d_all_lines[index*5+Y1] = points[j].y;
			d_all_lines[index*5+LEN] = length;

			index++;
		}
	}

	// Kept by root only, used to free up all points in min lines saved by root.
	min_line_points = (point_t*) allocate(sizeof(point_t) * my_line_count * 2);

	free(points);
}



void distrib_lines() {
	long l_base;
	int remainder;
	int *i_send_counts;
	int i_recv_doubs;
	int *i_displs;

	if(my_rank==ROOT) {
		// Base number of lines to send (lines being 5 doubles)
	 	l_base = my_line_count/nprocs;
	 	
	 	// if there are any remaining lines after the base amount is split up
	 	remainder = my_line_count%nprocs;

		// Amount of lines (5 doubles) to send to each process 
	 	i_send_counts = (int*) allocate(sizeof(int) * nprocs); 
	 	
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
	MPI_Scatterv(d_all_lines, i_send_counts, i_displs, MPI_DOUBLE, d_my_lines, 
		i_recv_doubs, MPI_DOUBLE, ROOT, MPI_COMM_WORLD);

	ln_my_lines = (line_t *) allocate(my_line_count*sizeof(line_t));

	double_array_to_struct(d_my_lines, ln_my_lines, my_line_count);

	if(my_rank == ROOT) {
		free(i_displs);
		free(i_send_counts);
		free(d_all_lines);
	}
}


/*
 * Builds the Greedy Triangulation. 
 */
void triangulate() {
	// The triangulation will be stored as an array of lines. The triangulation
	// is built on process zero iteratively as successive global minimal lines
	// are found. Allocate enough space to potentially hold every line.
	if(my_rank == ROOT) {
		triang = (line_t*) allocate(my_line_count*nprocs*sizeof(line_t));
	}

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

		int min_line_index = 0; // Will hold index of the global min line.

		if (my_unknown > 0) {
			// Convert this processes' minimal (smallest) line to an array of
			// five doubles for Allgather.
			double my_min_line[5];
			point_t p = *(ln_my_lines[0].p);   /// Make sure there is a line in ln_my_lines
			point_t q = *(ln_my_lines[0].q);
			
			my_min_line[X0] = p.x;
			my_min_line[Y0] = p.y;
			my_min_line[X1] = q.x;
			my_min_line[Y1] = q.y;
			my_min_line[LEN] = ln_my_lines[0].len;
			// Prepare an array to receive each processe's minimal line.
			double* recv_buf = (double*) allocate(5*nprocs*sizeof(double)); 
			// Make sure each process has an array of each processes' min line.
			MPI_Allgather(my_min_line, 5, MPI_DOUBLE, 
				          recv_buf, 5, MPI_DOUBLE, MPI_COMM_WORLD);

			// Find the global minimal line.
			for (int i = 0; i < nprocs; i++) {
				if (recv_buf[i*5+LEN] != -1) {	// If not an impossible line		
					// Compare the length of the current smallest line to the
					// i^th line's length. If the length is not positive ignore it
					// because it was a special value sent from a process with no
					// more lines of unknown status.
					if ((recv_buf[i*5+LEN] < recv_buf[min_line_index*5+LEN])) {
						min_line_index = i;
					}
					else if (recv_buf[min_line_index*5+LEN] == -1) {
						min_line_index = i;
					}
				} // end if not impossible
			} // end for

			// Will hold the minimal line.
			line_t* min_line = (line_t*) allocate(sizeof(line_t)); // // // // // // // // // // // // // // // // // // // // // // // // // // // // // // 
			
			// See if this process has the global min line, if so we must
			// adjust its number of lines of unknown status and set min_line.
			if (my_rank == min_line_index) {
				my_unknown--;
				min_line = &ln_my_lines[0];
			// Otherwise we must build the min_line from data in the recv_buf
			// (While making sure to use the appropriate index!)
			}
			else {
				// Get the minimal line
				point_t *p = (point_t*) allocate(sizeof(point_t));
				point_t *q = (point_t*) allocate(sizeof(point_t));
				p->x = recv_buf[min_line_index*5 + X0];
				p->y = recv_buf[min_line_index*5 + Y0];
				q->x = recv_buf[min_line_index*5 + X1];
				q->y = recv_buf[min_line_index*5 + Y1];

				min_line->p = p;
				min_line->q = q;
				min_line->len = recv_buf[min_line_index*5 + LEN];
				// free(p);																		// FREE THESE
				// free(q);																		// FREE THESE
			}

			// Have process zero add min_line to the triangulation.
			if (my_rank == ROOT) {
				triang[tlines] = *min_line;
				tlines++;
			}
			// Free the receive buffer
			free(recv_buf);
			// Allocate an array of lines to hold the lines that don't 
			// intersect with the global minimum.
			line_t* temp = (line_t*) allocate(my_line_count*sizeof(line_t));
			int end = my_unknown;
			int start = 0;
			if (my_rank == min_line_index) {
				end = my_unknown+1;
				start = 1;
			}
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
			// Have everyone but the root deallocate space associated with min_line
			if ((my_rank != ROOT) && (my_rank != min_line_index)) {
				// free(min_line->p);
				// free(min_line->q);
				free(min_line);
			}
		} // end if (my_unknown > 0)

		// If this process has no more lines of unknown status then it must still
	    // participate in global communications to avoid deadlock. Have it send
		// a special value to the other processes (which they will ignore). The 
		// special value is a line whose endpoints are at the origin and it has a
		// distance of -1.
		else {
			// Prepare an array to receive each processe's minimal line.
			double* recv_buf = (double*) allocate(5*nprocs*sizeof(double));
			MPI_Allgather(IMPOSSIBLE_LINE, 5, MPI_DOUBLE, 
				          recv_buf, 5, MPI_DOUBLE, MPI_COMM_WORLD);

			// If ROOT, then make the minimal line into a line struct and 
			// store it in the triagulation.
			if (my_rank == ROOT) {

				// Find the global minimal line.
				for (int i = 0; i < nprocs; i++) {
					if (recv_buf[i*5+LEN] != -1) {	// If not an impossible line		
						// Compare the length of the current smallest line to the
						// i^th line's length. If the length is not positive ignore it
						// because it was a special value sent from a process with no
						// more lines of unknown status.
						if ((recv_buf[i*5+LEN] < recv_buf[min_line_index*5+LEN])) {
							min_line_index = i;
						}
						else if (recv_buf[min_line_index*5+LEN] == -1) {
							min_line_index = i;
						}
					} // end if not impossible
				} // end for

				// Will hold the minimal line.
				line_t* min_line;

		
				// Get the minimal line
				min_line = (line_t*) allocate(sizeof(line_t)); 
				point_t *p = (point_t*) allocate(sizeof(point_t));
				point_t *q = (point_t*) allocate(sizeof(point_t));
				p->x = recv_buf[min_line_index*5 + X0];
				p->y = recv_buf[min_line_index*5 + Y0];
				q->x = recv_buf[min_line_index*5 + X1];
				q->y = recv_buf[min_line_index*5 + Y1];
				min_line_points[mlp_index] = *p;
				min_line->p = &min_line_points[mlp_index];
				mlp_index++;

				min_line_points[mlp_index] = *q;
				min_line->q = &min_line_points[mlp_index];
				mlp_index++;

				min_line->len = recv_buf[min_line_index*5 + LEN];


				// free(p);
				// free(q);
				// // Add line to triagulation
				triang[tlines] = *min_line;
				tlines++;
			}
			// Check if all of the lines have distance of -1. If so then we are done!
			int count = 0;
			for (int i = 0; i < nprocs; i++) {
				if (recv_buf[i*5 + LEN] == -1) {
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
		printf("Usage %s <input filename>\n", argv[0]);
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

	MPI_Barrier(MPI_COMM_WORLD);
	START_TIMER(generate)
	if (my_rank == ROOT) {
		gen_lines();
	}

	distrib_lines();
	MPI_Barrier(MPI_COMM_WORLD);
	STOP_TIMER(generate)

	  //                                      //
	 //  Sort the lines from small to large  //
	//                                      //
	MPI_Barrier(MPI_COMM_WORLD);
	START_TIMER(sort)
	qsort(ln_my_lines, my_line_count, sizeof(line_t), compare);
	MPI_Barrier(MPI_COMM_WORLD);
	STOP_TIMER(sort)

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
		fprintf(write_file, "%ld\n", (tlines-1));
		
		// Write triangulation to file.
		for (int i = 0; i < tlines-1; i++)
		{
		    point_t p = *(triang[i].p);
		    point_t q = *(triang[i].q);
		    fprintf(write_file, "(%lf, %lf) (%lf, %lf)\n", p.x, p.y, q.x, q.y);
		    // free(triang[i].p);
		    // free(triang[i].q);
		}
		
		fclose(write_file);
	}
	
	// Clean up and exit
	free(points_to_free);
	free(ln_my_lines);

	if (my_rank == ROOT) {
		free(triang);
		free(min_line_points);
	}

	MPI_Finalize();
	return (EXIT_SUCCESS);
}
