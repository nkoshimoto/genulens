CC = clang++
CFLAGS  = -g -O3
INCP = -lgsl

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
	$(CC) $(CFLAGS) -o genulens genulens.o option.o $(INCP)

# To create the object file option.o, we need the source
# files option.cpp and option.h:
#
option.o:  option.cpp option.h
	$(CC) $(CFLAGS) -c option.cpp

# To create the object file genulens.o, we need the source file
# genulens.cpp:
#
genulens.o:  genulens.cpp
	$(CC) $(CFLAGS) -c genulens.cpp $(INCP)

# To start over from scratch, type 'make clean'.  This
# removes the executable file, as well as old .o object
# files and *~ backup files:
#
clean: 
	$(RM) count *.o *~
