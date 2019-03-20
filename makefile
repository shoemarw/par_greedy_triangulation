default: wrappers.o greedy_triangle.o turtle.o generate_points.o par_greedy_triangle.o
	gcc --std=c99 -o tri greedy_triangle.o wrappers.o turtle.o -lm
	gcc --std=c99 -o gen generate_points.o wrappers.o turtle.o -lm
	mpicc --std=c99 -o par_tri par_greedy_triangle.o wrappers.o turtle.o -lm
	
generate_points.o: generate_points.c
	gcc -g --std=c99 -Wall -c generate_points.c
	
turtle.o: turtle.c
	gcc -g --std=c99 -Wall -c turtle.c

wrappers.o: wrappers.c
	gcc -g --std=c99 -Wall -c wrappers.c
	
greedy_triangle.o: greedy_triangle.c
	gcc -g  --std=c99 -Wall -c greedy_triangle.c

par_greedy_triangle.o:
	mpicc -g --std=c99 -Wall -c par_greedy_triangle.c
	
clean:
	rm *.o; rm tri; rm gen; rm par_tri
