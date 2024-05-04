# Define the C++ compiler
CXX = g++

# Compilation flags: -Wall for compiler warnings, -g for debugging information
CXXFLAGS = -Wall -g -Wno-sign-compare -Wno-unused-variable

# Include directories for header files
INCLUDES = -I/usr/include

# Library paths, if Armadillo is not in a standard location
LFLAGS = -L/usr/lib

# Libraries to link against, including Armadillo
LIBS = -larmadillo

# Source files
MO_POINTS_SRCS = generate_mo_points.cpp CNDO.cpp util.cpp AO.cpp

# Object files, derived from SRC3S
MO_POINTS_OBJS = $(MO_POINTS_SRCS:.cpp=.o)

# Executable names
MO_POINTS_MAIN = mo_points_main

# Targets
.PHONY: all depend clean

all:    $(MO_POINTS_MAIN)
	@echo  Orbital points executable named $(MO_POINTS_MAIN) has been compiled

$(MO_POINTS_MAIN): $(MO_POINTS_OBJS)
	$(CXX) $(CXXFLAGS) $(INCLUDES) -o $(MO_POINTS_MAIN) $(MO_POINTS_OBJS) $(LFLAGS) $(LIBS)

# Rule to convert .cpp to .o
.cpp.o:
	$(CXX) $(CXXFLAGS) $(INCLUDES) -c $<  -o $@

clean:
	$(RM) *.o *~ $(MO_POINTS_MAIN) *.txt

depend: $(MO_POINTS_SRCS)
	makedepend $(INCLUDES) $^

# DO NOT DELETE
