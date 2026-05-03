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
int split(char *splitw, const char *s, char *word[])
/* Count the number of words in a string and return indivual words */
{
   int count = 0;
   int  index;
   char wd[100];
   while (*s != '\0') {
      while (*s == *splitw)   /* skip white space */
           ++s;
      if (*s != '\0') {          /* found a word */
         ++count;
         index = 0;
         while (*s!=*splitw && *s != '\0') {
            wd[index] = *s;      /* extract the word */
            ++index;
            ++s;
         }
         wd[index] = '\0';    /* put word in list */
         //word[count-1] = new char[strlen(wd) + 1]; // for C++
         //word[count-1] = calloc(strlen(wd)+1, sizeof(char)); // for C
         word[count-1] = (char*)calloc(strlen(wd)+1, sizeof(char));
         strcpy(word[count-1], wd);
      }
   }
   return count;
}
