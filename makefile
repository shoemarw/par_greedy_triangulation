default: wrappers.o greedy_triangle.o turtle.o gen.o
	gcc --std=c99 -o triangulation greedy_triangle.o wrappers.o turtle.o -lm
	# gcc --std=c99 -o gen gen.o wrappers.o turtle.o -lm
	
gen.o: generate_points.c
	gcc -c generate_points.c
	
turtle.o: turtle.c
	gcc -c turtle.c

wrappers.o: wrappers.c
	gcc -c wrappers.c
	
greedy_triangle.o: greedy_triangle.c
	gcc -c greedy_triangle.c
	
clean:
	rm *.o; rm triangulation; rm gen
