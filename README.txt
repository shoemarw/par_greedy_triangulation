GNU General Public License v3.0
//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                         USING THE SOFTWARE                           //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

To use the software first use the makefile, three executables are made.
Before using the makefile make sure to load mpi. The makefile produces the
following executables: tri, gen, and par_tri. tri triangluates a point set
serially and par_tri triangulates a point set using MPI. The executables are
described below, the other files are too.

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                          THE EXECUTABLES                             //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

gen:
gen generates an output file containing points. The output file is
made to be used as input to both of tri and par_tri. This executable takes 
two command line parameters: <num points> and <seed>. <num points> is the
number of points that you wish to generate. The points are generated by
a pseudo random number generator (a linear congruential generator (LCG)),
thus the command line argument <seed> is the seed for the LCG. The output
file has the following format: Line one is one number which specifies the
number of points of in the input file. Each subsequent line contains two
numbers separated by a space; each line in the file corresponds to a point.

tri:
tri takes an input file containing a point set; the format of the input
file is specified in the description of the gen executable. tri produces an 
output file containing the the greedy triangulation of the point set
specified in its input file. It may also produce a .bmp image of the 
triangulation if IMAGE is defined in greedy_triangulation.c when the 
makefile is used. The output file (not the .bmp file) containing the 
triangulation has the following format: the first line contains a single
number representing the number of lines in the triangulation. Each subsequent
line contains two tuples representing the endpoints of a line making the 
triangulation.

par_tri:
same as tri except this is the parallel version. Uses MPI.

//////////////////////////////////////////////////////////////////////////
//                                                                      //
//                        DESCRIPTION OF FILES                          //
//                                                                      //
//////////////////////////////////////////////////////////////////////////

generate_points.c:
This code in this file is used to create input files for the triangulation
executables. Its functionality is described where the gen executable is
described above.

greedy_triangle.c:
This code is the serial version of the greedy triangulation algorithm. Its
functionality is described where the tri executable is described above.

greedy_triangle.h:
This header file contains declarations of structures, helper functions and 
wrappers used in greedy_triangle.c, generate_points.c, and
par_greedy_triangle.c

makefile:
Described above in the "USING THE SOFTWARE" section.

par_greedy_triangle.c:
Same as greedy_triangle.c except parallelized by using MPI.

test20pts31:
This is a sample input file for either the tri or the par_tri executables.
This is also a sample output of generate_points.c (the gen executable).
to use this with try type:
    ./tri test20pts31
This was produced by gen using the command line parameters: 20 for number
of points and 31 as the seed for the LCG.

turtle.c:
Takes care of creating images: produced by Dr. Lam.

turtle.h:
Header file for turtle.c.

wrappers.c:
Contains the definitions for helper functions and wrappers.
