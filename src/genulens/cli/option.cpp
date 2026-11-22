// This program was originally written by Ian A. Bond and modified by Takahiro Sumi afterwards. The copyright belong to them.
#include "option.h"

//------------------------------------------------------------------------
char* getOptions(int argc, char *argv[], const char *argname, int argno, char *def_word = (char *)"undef")
{
  int i;
  for (i=1;i<argc;i++){
     if (!strcmp(argv[i], argname)) break;
  }
  //printf ("%d\n", strlen(argv[i+argno]));
  //printf ("%d %d %d\n",argc, i, argno);
  if (i+argno<argc) return argv[i+argno];
  else return def_word;
}
//------------------------------------------------------------------------
int getOptioni(int argc, char *argv[], const char *argname, int argno, int def_value)
{
   char *word = getOptions(argc, argv, argname, argno, (char *)"undef");
   if (!strcmp(word, "undef")) return def_value;
   else return atoi(word);
}
//------------------------------------------------------------------------
double getOptiond(int argc, char *argv[], const char *argname, int argno, double def_value)
{
   char *word = getOptions(argc, argv, argname, argno, (char *)"undef");
   if (!strcmp(word, "undef")) return def_value;
   else return atof(word);
}
/* -------------------------------------------------------------------*/
void split_table_fields(char *line, std::vector<const char *> &fields)
{
    fields.clear();
    while (*line) {
        line += strspn(line, " ");
        if (!*line) break;
        fields.push_back(line);
        char *end = strchr(line, ' ');
        if (!end) break;
        *end = '\0';
        line = end + 1;
    }
}
