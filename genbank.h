#ifndef _libgenbank
#define _libgenbank


#define MAXCOL 1000

typedef struct {
    char **str;     //the PChar of string array
    size_t num;     //the number of string
}IString;

char buf[MAXCOL];
char line[MAXCOL];
IString subline;
IString data;
IString location;

int split(char *src, char *delim, IString* istr);
void get_CDS(char name[]);

#endif //_libgenbank
