/*
 * greedy_triangle.c
 * 
 * Serial version of the greedy triangulation algorithm.
 * We read a file containing a list of points. The input file must
 * have a number on its first line indicating the number of lines
 * in the file. Each line there after must contain two numbers 
 * separated by a space. Each line is the Cartesian representation
 * of a point in in the standard x-y plane.
 *
 * CS 470 Research Project.
 * Original serial version.
 *
 * Name(s): Randy, Eliza, Alex
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
#define IMAGE

int main(int argc, char *argv[]) {
	
	// Make sure we get the expected input.
	if (argc != 2) {
		printf("Usage %s <filename>, argv[0] \n", argv[0]);
		exit(EXIT_FAILURE);
	}
	
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
	
	  //                      //
	 //  Generate all lines  //
	//                      //
	
	// utility struct for timing calls
    struct timeval tv;
	START_TIMER(generate)
	// Make all possible line segments between the points
	// and compute the length of each line.
	// Compute the number of lines, note that the below formula
	// will always resolve as an int because one of num_points or
	// num_points-1 will be divisible by 2.
	int num_lines = ((num_points)*(num_points-1))/2;
	line_t* lines = (line_t*) allocate(num_lines * sizeof(line_t));
	
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
			lines[index] = *l;
			index++;
			free(l);
		}
	}
	STOP_TIMER(generate)
	
	  //                                      //
	 //  Sort the lines from small to large  //
	//                                      //
	
	START_TIMER(sort)
	qsort(lines, num_lines, sizeof(line_t), compare);
	STOP_TIMER(sort)
		
      //                                   //
     //  Greedily build the tringulation  //
    //	                                 //
	
	START_TIMER(triangulate)
	// The triangulation will be stored as an array of lines.
	// Allocate space for the triangulation.
	line_t* triang = (line_t*) allocate(num_lines*sizeof(line_t));
	
	// unknown keeps track of how many lines remain that may or may not be in the 
	// triangulation. As soon as a line is found to be in the triangulation it is
    // is decremented. As soon as it is known that a particular line does not 
	// belong to the triangulation it is decremented.
	long unknown = num_lines;
	long tlines  = 0;        // Tracks the number of lines in the triangulation.
	// Keep going until there are no more lines of unknown status.
	while (unknown > 0) {
		// Add the smallest element of lines array to the triangulation.
		triang[tlines] = lines[0]; // lines[0] is always the smallest line.
		tlines++; // increment the count of lines in triang.
		unknown--; // decrement unknown since a line has been added to triang.
		line_t* temp = allocate(num_lines * sizeof(line_t));
		// Add all remaining lines to temp that don't intersect with the newly
        // added line to triang. (this is the step where we 'get rid' of lines 
		// which conflict with the line newly added to the triagulation).
		int end = unknown;
		int temp_size = 0;
		
		for (int j = 1; j < end + 1; j++)
		{
		    // Run intersection test. Note triang[tlines-1] is the line just added
		    // to the triangulation. If lines[i] does not intersect the newest 
		    // member of triang then it could be part of the triangulation. If 
		    // it does intersect then don't add it to temp because it cant be part
		    // of the triangulation. If two lines share an endpoint then we consider
		    // them to NOT intersect for the sake of the triangulation.
		    if (share_endpoint(&lines[0], &lines[j]) ||
		         !intersects(&lines[0], &lines[j]))
		    {
		        temp[temp_size] = lines[j];
		        temp_size++;
		    }
		    else
		    {
		        unknown--;
		    }
		}
		
		// Write all of the valid lines from temp to lines[] so we can iterate again.		
		copy_array(temp, lines, temp_size);
		free(temp);
	}
	STOP_TIMER(triangulate)
	
	  //                                     //
     //  Triangulation Done, Display Stats  //
    //	                                   //
	
	// These stats are only for the portions of the code specific to the three
	// phases of building the greedy triangulation. Generate all lines, sort the
	// lines in non-decreasing order, and greedily adding line segments to the
	// triangulation.
	printf("Gent: %.4f  Sort: %.4f  Tria: %.4f\n",
	        GET_TIMER(generate), GET_TIMER(sort), GET_TIMER(triangulate));
	
	// Generate image
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
	
	// Clean up and exit
	free(triang);
	free(points);
	free(lines);
	return (EXIT_SUCCESS);
}
