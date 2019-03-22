/*
 * greedy_triangle.h
 *
 * Header file for declaration of structures, helper 
 * functions, and wrapper functions
 * 
 * CS 470 Research Project.
 * Original serial version.
 *
 * Name(s): Randy, Eliza, Alex
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>

#include "turtle.h"

/*
 * Represents a point in the plane using cartesian coordinates.
 */
typedef struct
{
    double x;
    double y;
} point_t;

/*
 * Represents a line segment between two points.
 */
typedef struct
{
    point_t *p;
    point_t *q;
    double len;
} line_t;

void print_point(point_t* point);

void print_line(line_t* line);

bool is_equal(point_t* p, point_t* q);

double distance(point_t* p, point_t* q);

int compare(const void* a, const void* b);

double max(double m, double n);

double min(double m, double n);

bool lies_on(point_t* p, point_t* q, point_t* r);

int orient(point_t* p, point_t* q, point_t* r);

int intersects(line_t* a, line_t* b);

bool share_endpoint(line_t* a, line_t* b);

void copy_array(line_t from[], line_t to[], int size);

void* allocate(size_t size);

FILE* open_file(char* filename, char* operation);

void error(char* message);

void generate_image(line_t lines[], int number_lines);
