SOURCES =	main.cpp MatrixClass.cpp UnSparseClass.cpp SparseClass.cpp 

OBJECTS = $(SOURCES: .cpp = .o)

EXECUTABLE = Matrix

SHELL = /bin/sh

CC := gcc
CXX := g++ 
FC = gfortran
#CFLAGS := -O0 -g -std=c99 -Wall -Wextra -Wshadow -pedantic -Werror
CXXFLAGS := -O0 -g -std=c++11 -Wall -Wextra -Wshadow -pedantic -Weffc++ -Werror
FFLAGS = -Wall -g

APPLECFLAGS = -m64 -arch x86_64

FRAMEWORK = #-F/System/Library/Frameworks/Accelerate.framework -framework Accelerate

APPLELFLAGS =  -L/opt/local/lib -llapack -L/opt/local/lib -larpack \
	/usr/local/lib/x86_64/libgfortran.3.dylib 

#LIBDIR = /opt/intel/Compiler/11.1/038/mkl/lib/em64t

#MKLFLAGS = -L${LIBDIR} -lmkl_lapack 

### For Mac OS X
LIBS = 	$(FRAMEWORK) $(APPLELFLAGS)
#CXX = g++

### For Nautilus
#LIBS = $(ARPACK_LIB) $(GSL_LIB) $(MKL_LIB) $(MKLFLAGS) -llapack
         

all : 	$(SOURCES) $(EXECUTABLE) 
	$(CXX) $(CXXFLAGS) $(LIBS) $(OBJECTS) -o $(EXECUTABLE)

#GeneralFunctions.o : GeneralFunctions.cpp
#	$(CXX) -static -c GeneralFunctions.cpp -lgsl -lgslcblas -lm

$(EXECUTABLE) :	$(OBJECTS)
	$(CXX) $(CXXFLAGS) $(LIBS) $(OBJECTS) -o $(EXECUTABLE)

.cpp.o:
	$(CXX) $(CXXFLAGS) $< -o $(EXECUTABLE)

.f.o:
	$(FC) $(FFLAGS) $< -o $(EXECUTABLE)

clean: 
	rm -f *.o *~ $(EXECUTABLE) 




