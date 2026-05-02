# This is an example of Makefile.
# You might be able to use this Makefile without any edition.
# If the default form does not work, the C++ compiler specified in the CC line or the GSL paths specified in the INCLUDE and LINK lines are probably inappropriate.
# Please replace them with what you have.
# The paths for GSL can be found by 
#   $ gsl-config --libs   (<- for LINK)
#   $ gsl-config --cflags (<- for INCLUDE)
# If the command gsl-config does not work in the terminal, it probably means that the GSL lib is not installed, or unknown to the OS.
#
CC = clang++
#CC = g++-mp-9
CFLAGS  = -g -O3
# CFLAGS  = -g
LIBS = -lm -lgsl -lgslcblas
INCLUDE = -I/opt/local/include
LINK = -L/opt/local/lib

# typing 'make' will invoke the first target entry in the file 
# (in this case the default target entry)
# you can name this target entry anything, but "default" or "all"
# are the most commonly used names by convention
#
default: genulens

# To create the executable file genulens we need the object files
# genulens.o and option.o:
#
genulens: genulens.o option.o
	$(CC) $(CFLAGS) -o genulens genulens.o option.o $(LINK) $(LIBS)

# To create the object file option.o, we need the source
# files option.cpp and option.h:
#
option.o:  option.cpp option.h
	$(CC) $(CFLAGS) -c option.cpp

# To create the object file genulens.o, we need the source file
# genulens.cpp:
#
genulens.o:  genulens.cpp
	$(CC) $(CFLAGS) -c genulens.cpp $(INCLUDE)

# To start over from scratch, type 'make clean'.  This
# removes the executable file, as well as old .o object
# files and *~ backup files:
#
clean: 
	$(RM) count *.o *~
