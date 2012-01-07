.PHONY: default clean

CC = g++

default: spline.o


#spline ========================================
spline.o: spline.C
spline.o: spline.H
	${CC} ${CCFLAGS} -c -o spline.o spline.C


clean:
	rm -rf  *.o *~
