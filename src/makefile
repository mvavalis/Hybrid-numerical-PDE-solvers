.PHONY:  default opt debug  clean clean_output

CC = g++

moredir = more

splinedir = Burkardt/spline
splineobj = ${splinedir}/spline.o

MBAdir = mba-1.1
MBAlibdir = ${MBAdir}/lib
MBAinclude = ${MBAdir}/include
MBAo = MBA.o MBAdata.o UCBsplineSurface.o UCBsplines.o UCButils.o
MBAobj = $(patsubst %, ${MBAlibdir}/%, ${MBAo}) 

DEALIIdir = /home/gesareli/wd/deal.II
#DEALIIdir = /home/manolis/pdd/deal.II
DEALIIinclude = $(patsubst %, ${DEALIIdir}/%, deal.II/include base/include lac/include contrib/boost/include)
DEALIIlibdir = ${DEALIIdir}/lib
DEALIIlibfiles = libbase  libdeal_II_1d  libdeal_II_2d  libdeal_II_3d  liblac

LINK = -lpthread -lm `pkg-config --libs --cflags gsl`

default: DEALIIlib = $(patsubst %, ${DEALIIlibdir}/%.so, ${DEALIIlibfiles}) 
default: CCFLAGS = -Wall
#default: CCFLAGS = -Wall -Wno-unused-variable
default: MBA_MODE = nopt
default: main.o

opt: DEALIIlib = $(patsubst %, ${DEALIIlibdir}/%.so, ${DEALIIlibfiles}) 
opt: CCFLAGS = -O2 -Wall
#opt: CCFLAGS = -O2 -Wall -Wno-unused-variable
opt: MBA_MODE = opt
opt: main.o

debug: DEALIIlib = $(patsubst %, ${DEALIIlibdir}/%.g.so, ${DEALIIlibfiles}) 
debug: CCFLAGS = -pg -Wall
debug: MBA_MODE = nopt
debug: main.o


#main ===========================================
main.o: main.cpp ${splineo}
main.o: dd.h Problem.h Pdd.h MCDriver.h Inter3DDriver.h LaplaceDriver.h LaplaceSolve.h ${moredir}/Ranq1.h ${moredir}/Vdcbin.h
	@make lib MODE=${MBA_MODE} --directory=${MBAdir} -f Makefile CCFLAGS=${CCFLAGS}
	@make --directory=${splinedir} -f makefile CCFLAGS=${CCFLAGS}
	${CC} ${CCFLAGS} -o main.out main.cpp $(patsubst %, -I%, ${DEALIIinclude}) -I ${MBAinclude} ${splineobj} ${MBAobj} ${DEALIIlib} ${LINK}

#more ==========================================
clean:
	@make --directory=${MBAdir} -f Makefile clean
	@make --directory=${splinedir} -f makefile clean
	rm -rf *.o *~ main.out gmon.out

clean_output:
	rm -rf ../solution/*.gnuplot ../solution/*.tecplot ../solution/*.dx ../solution/*.povray

