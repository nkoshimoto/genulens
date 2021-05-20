CC = clang++
CFLAGS  = -g -O3

# typing 'make' will invoke the first target entry in the file 
# (in this case the default target entry)
# you can name this target entry anything, but "default" or "all"
# are the most commonly used names by convention
#
default: genulens

# To create the executable file count we need the object files
# countwords.o, counter.o, and scanner.o:
#
genulens: genulens.o option.o
	$(CC) $(CFLAGS) -o genulens genulens.o option.o

# To create the object file countwords.o, we need the source
# files countwords.c, scanner.h, and counter.h:
#
option.o:  option.cpp option.h
	$(CC) $(CFLAGS) -c option.cpp

# To create the object file counter.o, we need the source files
# counter.c and counter.h:
#
genulens.o:  genulens.cpp
	$(CC) $(CFLAGS) -c genulens.cpp

# To start over from scratch, type 'make clean'.  This
# removes the executable file, as well as old .o object
# files and *~ backup files:
#
clean: 
	$(RM) count *.o *~
