/*
 * wrappers.c
 *
 * Declaration of helper functions and wrappers.
 * 
 * CS 470 Research Project.
 * Original serial version.
 *
 * Name(s): Randy, Eliza, Alex
 */
 
#include "greedy_triangle.h"


/*
 * Prints a point.
 */
void print_point(point_t* point)
{
    printf(" (%f , %f) ", point->x, point->y);
}

/*
 * Prints a line.
 */
void print_line(line_t* line)
{

    print_point(line->p);
	print_point(line->q);
	printf("\n");
}

/*
 * Determines if two points are equal. Retruns true if 
 * they are, false otherwise.
 */
bool is_equal(point_t* p, point_t* q)
{
    return p->x == q->x && p->y == q->y;
}

/* 
 * Computes the Euclidian distance between two points.
 */
double distance(point_t* p, point_t* q)
{
    double delta_x = p->x - q->x;
	double delta_y = p->y - q->y;
	return sqrt(delta_x * delta_x + delta_y * delta_y);
}

/* 
 * Compares two lines. Used for sorting lines with qsort.
 */
int compare(const void* a, const void* b)
{
    double result = ((line_t*) a)->len - ((line_t*) b)->len;
	// Make sure that cmp returns an int.
	if (result < 0) {
		return -1;
	}
	else if (result > 0) {
		return 1;
	}
	else {
		return 0;
	}
}

/*
 *  Returns the larger of m,n.
 */
double max(double m, double n)
{
    if (m > n) {
		return m;
	}
	else {
		return n;
	}
}

/*
 * Returns the smaller of m,n.
 */
double min(double m, double n)
{
    if(m < n) {
		return m;
	}
	else {
		return n;
	}
}

/*
 * Helper function for intersects():
 * Only use this function for colinear points p,q,r. Returns true (1)
 * if q lies on the line segment formed by p&r. adapted from:
 * https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
 */
bool lies_on(point_t* p, point_t* q, point_t* r)
{
    return (q->x <= max(p->x, r->x) && 
	    q->x >= min(p->x, r->x) &&
		q->y <= max(p->y, r->y) &&
		q->y >= min(p->y, r->y));
}

/*
 * Helper function for intersects():
 * Suppose p, q, and r form the verticies of a triangle. This
 * function informs us of the following: If  we traverse the
 * triangle starting at p then going to q then going to r and
 * finally back to r; this function tells us if the traversal
 * was clockwise, counter-clockwise, or if the points are 
 * colinear. 0 = colinear. 1 = clockwise. 2 = counter-clockwise.
 */
int orient(point_t* p, point_t* q, point_t* r)
{
    int result = (q->y - p->y) * (r->x - q->x) -
	             (q->x - p->x) * (r->y - q->y);
    
	// Check for colinear.
    if (result == 0) {
		return 0; 
	}
	// Check for clockwise
	else if (result > 0) {
		return 1;
	}
	// Otherwise its counter-clockwise.
    else {
		return 2;
	}
}

/*
 * Returns true (1) if the two line segments pointed to by r
 * and s intersect. Returns false (0) otherwise. The logic of
 * this function was adapted from code found at the url:
 * https://www.geeksforgeeks.org/check-if-two-given-line-segments-intersect/
 */
int intersects(line_t* a, line_t* b)
{
    // Get the points associated with the lines.
	point_t *f_pt1 = a->p;
	point_t *f_pt2 = a->q;
	point_t *g_pt1 = b->p;
	point_t *g_pt2 = b->q;
	// Find the orientations of all possible subsets with 3 points.
    int orient1 = orient(f_pt1, f_pt2, g_pt1); 
    int orient2 = orient(f_pt1, f_pt2, g_pt2); 
    int orient3 = orient(g_pt1, g_pt2, f_pt1); 
    int orient4 = orient(g_pt1, g_pt2, f_pt2); 
  
    // See if an intersection occured. 
    if (orient1 != orient2 && orient3 != orient4) 
        return 1; 
	
      //                                                    //
     // Cases where some subset of 3 points are colinear.  //
	//                                                    // 
	
    // See if f_pt1, f_pt2 and g_pt1 are colinear and f_pt2 lies on
	// the line between f_pt1 and g_pt1.
    if (orient1 == 0 && lies_on(f_pt1, g_pt1, f_pt2)){
		return 1; 
	}
  
    // See if f_pt1, f_pt2 and g_pt1 are colinear and f_pt2 lies on
	// the line between f_pt1 and g_pt2.
    if (orient2 == 0 && lies_on(f_pt1, g_pt2, f_pt2)) {
		return 1; 
	}
  
    // See if g_pt1, g_pt1 and f_pt1 are colinear and g_pt2 lies on 
	// the line between g_pt1 and f_pt1. 
    if (orient3 == 0 && lies_on(g_pt1, f_pt1, g_pt2)) {
		return 1;
	}		
  
    // See if g_pt1, g_pt1 and f_pt2 are colinear and g_pt2 lies on
    // the line between g_pt1 and f_pt2.	 
    if (orient4 == 0 && lies_on(g_pt1, f_pt2, g_pt2)) {
		return 1;
	}		
  
    // No intersection detected.
    return 0;
}

/*
 * Determines if two lines share an endpoint.
 */
bool share_endpoint(line_t* a, line_t* b)
{
    return is_equal(a->p,b->p) || is_equal(a->p,b->q) ||
           is_equal(a->q,b->p) || is_equal(a->q,b->q);
}

/*
 * Copies 'size' values of 'from' to array 'to'
 */
void copy_array(line_t from[], line_t to[], int size)
{
    for (int i = 0; i < size; i++)
    {
        to[i] = from[i];
    }
}

/*
 * Wrapper function for calloc. It checks for errors as well.
 */
void* allocate(size_t size)
{
    void* address = malloc(size);
    if (!address)
    {
        fprintf(stderr, "Cannot malloc, out of memory\n");
        exit(EXIT_FAILURE);
    }
    
    memset(address, 0, size);
    return address;
}

/*
 * Opens a file with the given operation, also checks for error.
 */
FILE* open_file(char* filename, char* operation)
{
    FILE* file = fopen(filename, operation);
    if (!file)
    {
        fprintf(stderr, "ERROR: Could not open %s\n", filename);
        exit(EXIT_FAILURE);
    }
    return file;
}

/*
 * Prints to the standard error and terminates the program with failure.
 */
void error(char* message)
{
    fprintf(stderr, "%s", message);
    exit(EXIT_FAILURE);
}

/*
 * Generates a .bmp file image with the triangulation lines.
 */
void generate_image(line_t lines[], int number_lines)
{
    turtle_init(1920, 1080); // standard resolution (width, height)
    
    turtle_draw_line(-960, 0, 960, 0); // draws x-axis
    turtle_draw_line(0, -540, 0, 540); // draws y-axis
    
    turtle_set_pen_color(255, 0, 0); // sets the color of the pen to red
                                  // triangles will appear in red axis in black
    
    int scale = 10; // Scales all of the points to fit the window
                                  
    for(int i = 0; i < number_lines; i++)
    {
        point_t p = *(lines[i].p);
        point_t q = *(lines[i].q);
        
        int p_x = (int) p.x;
        int p_y = (int) p.y;
        
        int q_x = (int) q.x;
        int q_y = (int) q.y;
        
        turtle_draw_line(p_x * scale, p_y * scale , q_x * scale , q_y * scale);
    }
    
    turtle_save_bmp("triangulation_graphics.bmp"); // save image
    
    turtle_cleanup(); // clean turtle
    
}


