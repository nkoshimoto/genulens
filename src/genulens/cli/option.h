#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>
#include <vector>

char* getOptions(int argc, char *argv[], const char *argname, int argno, char *def_word);
int getOptioni(int argc, char *argv[], const char *argname, int argno, int def_value);
double getOptiond(int argc, char *argv[], const char *argname, int argno, double def_value);

// Fields point into line and are valid until the caller reuses that buffer.
void split_table_fields(char *line, std::vector<const char *> &fields);
