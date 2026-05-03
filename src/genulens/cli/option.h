#include <stdio.h>
#include <stdlib.h>
#include <ctype.h>
#include <string.h>

char* getOptions(int argc, char *argv[], const char *argname, int argno, char *def_word);
int getOptioni(int argc, char *argv[], const char *argname, int argno, int def_value);
double getOptiond(int argc, char *argv[], const char *argname, int argno, double def_value);

int split(char *splitw, const char *s, char *word[]);
