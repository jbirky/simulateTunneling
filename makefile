CC=g++

# LFLAGS = -lm /usr/share/doc/blas -lcblas

LFLAGS = -lm -L/System/Library/Frameworks/Python.framework/Versions/2.6/Extras/lib/python/scipy/lib/blas/ -lcblas 

simulateTunneling: simulateTunneling.cpp
	$(CC) -o simulateTunneling simulateTunneling.cpp $(LFLAGS); ./simulateTunneling
	