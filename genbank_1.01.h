#ifndef _libgenbank
#define _libgenbank

#include<stdlib.h>
#include<string.h>
#include<stdio.h>

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

int split(char *src, char *delim, IString* istr)//split src，将分割后的子串存入结构体的str中
{
    int i;
    char *str = NULL, *p = NULL;
 
    (*istr).num = 1;
	str = (char*)calloc(strlen(src)+1,sizeof(char));
	if (str == NULL) return 0;
    (*istr).str = (char**)calloc(1,sizeof(char *));
    if ((*istr).str == NULL) return 0;
    strcpy(str,src);
 
	p = strtok(str, delim);
	(*istr).str[0] = (char*)calloc(strlen(p)+1,sizeof(char));
	if ((*istr).str[0] == NULL) return 0;
 	strcpy((*istr).str[0],p);
	for(i=1; p = strtok(NULL, delim); i++)
    {
        (*istr).num++;
        (*istr).str = (char**)realloc((*istr).str,(i+1)*sizeof(char *));
        if ((*istr).str == NULL) return 0;
        (*istr).str[i] = (char*)calloc(strlen(p)+1,sizeof(char));
        if ((*istr).str[0] == NULL) return 0;
        strcpy((*istr).str[i],p);
    }
    free(str);
    str = p = NULL;
 
    return 1;
}

void get_CDS(char name[])
{
    FILE *fp = fopen(name,"r");
    while(fgets(line,MAXCOL,fp))
    {
        split(line," ",&subline);
        if (strcmp(subline.str[0],"CDS") == 0)
            break;
        else free(subline.str);
    }

    split(subline.str[1],"<)",&data);
    free(subline.str);
    split(data.str[1],".",&location);
    free(data.str);

    int begin = atoi(location.str[0]);
    int end = atoi(location.str[1]);
    free(location.str);

    int line_begin = (begin - 1) / 60;
    int line_end = (end - 1) / 60;
    int col_begin = (begin - 1) % 60;
    int col_end = (end - 1) % 60;

    while(fgets(line,MAXCOL,fp))
    {
        split(line," ",&subline);
        if (strcmp(subline.str[0],"1") == 0)
            break;
        else free(subline.str);
    }
    
    int i = 0;
    do  //序列读取
    {
        if (i == line_begin){
            split(line," ",&subline); //读到最后
            int subcol = col_begin / 10 + 1; //获取字串所在的subline的位置
            int col_start = col_begin % 10;  //获取列数
            for (;col_start < 10;++col_start)
            {
                printf("%c",subline.str[subcol][col_start]);
            }
            subcol++;
            for (;subcol < 7;++subcol)
            {
                printf("%s",subline.str[subcol]);
            }
            free(subline.str);
        }

        else if (i > line_begin && i < line_end){ //输出整行
            int subcol = 1;
            split(line," ",&subline); //读到最后
            for (;subcol < 7;++subcol)
            {
                printf("%s",subline.str[subcol]);
            }
            free(subline.str);
        }

        else if (i == line_end){ //输出到col_end
            split(line," ",&subline); 
            int subcol = col_end / 10 + 1;
            int col_over = col_end % 10;
            
            int col = 1;
            for (;col < subcol;++col)
            {
                printf("%s",subline.str[col]);
            }
            int i = 0;
            for (;i <= col_over;++i)
            {
                printf("%c",subline.str[subcol][i]);
            }
            free(subline.str);
            break;
        }
        
        i++;
    }while(fgets(line,MAXCOL,fp));
    
}

/*void getline_n(FILE *Target,int n,char buf[])  //获取文件的第n行(fgets的第n - 1行)
{
    for (int i = 0;i < n - 1;++i)  //跳过前面的n - 1行
    {
        fscanf(Target,"%*[^\n]%*c");
    }

    fgets(buf,MAXCOL,Target);
}*/
#endif //_libgenbank
