default: wrappers.o greedy_triangle.o turtle.o gen
	gcc --std=c99 -o triangulation greedy_triangle.o wrappers.o turtle.o -lm
	
turtle.o: turtle.c
	gcc -c turtle.c

wrappers.o: wrappers.c
	gcc -c wrappers.c
	
greedy_triangle.o: greedy_triangle.c
	gcc -c greedy_triangle.c
	
gen: generate_points.c
	gcc --std=c99 -o generate_points.c wrappers.o turtle.o -lm
	
clean:
	rm *.o; rm triangulation
