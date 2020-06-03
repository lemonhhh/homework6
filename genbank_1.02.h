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
char *pro_euka_sequence;

int split(char *src, char *delim, IString* istr);
void Reverse(char s[]);
void tackle_euka(FILE *fp);
void tackle_proka(FILE *fp);
void get_sequence(int begin,int end,FILE *fp,int* i);
void get_sequence_p(int begin,int end,FILE *fp);
void reverse_basic(char group[],char basic,int *j);
void get_CDS(char name[]);

/*void getline_n(FILE *Target,int n,char buf[])  //获取文件的第n行(fgets的第n - 1行)
{
    for (int i = 0;i < n - 1;++i)  //跳过前面的n - 1行
    {
        fscanf(Target,"%*[^\n]%*c");
    }

    fgets(buf,MAXCOL,Target);
}*/

void get_CDS(char name[])
{
    FILE *fp = fopen(name,"r");
    while(fgets(line,MAXCOL,fp))   //读取文件直至CDS那一行
    {
        split(line," ",&subline);
        if (strcmp(subline.str[0],"CDS") == 0)
            break;
        else free(subline.str);
    }
    
    if (subline.str[1][0] == 'c')
        tackle_proka(fp);
    else if (subline.str[1][0] == 'j')
        tackle_euka(fp);
    
}
void tackle_proka(FILE *fp)    //处理原核的基因序列
{
    split(subline.str[1],"<)",&data);
    free(subline.str);
    split(data.str[1],".",&location);
    free(data.str);

    int begin = atoi(location.str[0]); //获取到相应的起始位点和终止位点
    int end = atoi(location.str[1]);
    free(location.str);

    char sequence[] = "1";
    while(fgets(line,MAXCOL,fp)) //定位到序列行 
    {
        split(line," ",&subline);
        if (strcmp(subline.str[0],sequence) == 0)
            break;
        else free(subline.str);
    }
    get_sequence_p(begin,end,fp);
    Reverse(pro_euka_sequence);
    printf(pro_euka_sequence);
}
void tackle_euka(FILE *fp)      //处理真核的基因序列
{
    split(subline.str[1],"()",&data);
    free(subline.str);
    split(data.str[1],".,",&location);
    free(data.str);

    int all_region = location.num;
    int *section_start = (int *)malloc(sizeof(int) * all_region/2);
    int *section_end = (int *)malloc(sizeof(int) * all_region/2);
    int j = 0;
    for (int i = 0; i < all_region; i +=2)  //把起始终止的片段存储在数组中
    {
        section_start[j] = atoi(location.str[i]);
        section_end[j++] = atoi(location.str[i + 1]);
    }
    free(location.str);

    char sequence[] = "1";
    while(fgets(line,MAXCOL,fp))  //将文件光标转移到序列起始的那一行
    {
        split(line," ",&subline);
        if (strcmp(subline.str[0],sequence) == 0)
            break;
        else free(subline.str);
    }
    int line_num_ = 0;
    for (int k = 0; k < all_region / 2;++k)
    {
        printf("Section %d:\n",k);
        get_sequence(section_start[k],section_end[k],fp,&line_num_);
    }
}

void get_sequence(int begin,int end,FILE *fp,int* i)    //处理真核生物的序列提取
{
    int line_begin = (begin - 1) / 60;
    int line_end = (end - 1) / 60;
    int col_begin = (begin - 1) % 60;
    int col_end = (end - 1) % 60;
    do  //序列读取
    {
        if ((*i) == line_begin){
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

        else if ((*i) > line_begin && (*i) < line_end){ //输出整行
            int subcol = 1;
            split(line," ",&subline); //读到最后
            for (;subcol < 7;++subcol)
            {
                printf("%s",subline.str[subcol]);
            }
            free(subline.str);
        }

        else if ((*i) == line_end){ //输出到col_end
            split(line," ",&subline); 
            int subcol = col_end / 10 + 1;
            int col_over = col_end % 10;
            
            int col = 1;
            for (;col < subcol;++col)
            {
                printf("%s",subline.str[col]);
            }
            int k = 0;
            for (;k <= col_over;++k)
            {
                printf("%c",subline.str[subcol][k]);
            }
            free(subline.str);
            break;
        }
        
        (*i)++;
    }while(fgets(line,MAXCOL,fp));
    printf("\n");
}
void get_sequence_p(int begin,int end,FILE *fp) //原核生物基因 CDS序列输出
{
    int line_begin = (begin - 1) / 60;
    int line_end = (end - 1) / 60;
    int col_begin = (begin - 1) % 60;
    int col_end = (end - 1) % 60;
    pro_euka_sequence = (char *)malloc(sizeof(char) * (end - begin + 300)); //将序列存储到字符串中(还有换行符)
    int realnum = 0;
    int i = 0;
    int z; //计数因子
    do  //序列读取
    {
        if (i == line_begin){
            split(line," ",&subline); //读到最后
            int subcol = col_begin / 10 + 1; //获取字串所在的subline的位置
            int col_start = col_begin % 10;  //获取列数
            for (;col_start < 10;++col_start)
            {
                reverse_basic(pro_euka_sequence,subline.str[subcol][col_start],&realnum);
            }
            subcol++;
            for (;subcol < 7;++subcol)
            {
                for (z = 0;z < 10;++z)
                    reverse_basic(pro_euka_sequence,subline.str[subcol][z],&realnum);
            }
            pro_euka_sequence[realnum++] = '\n';
            free(subline.str);
        }

        else if (i > line_begin && i < line_end){ //输出整行
            int subcol = 1;
            split(line," ",&subline); //读到最后
            for (;subcol < 7;++subcol)
            {
                for (z = 0;z < 10;++z)
                    reverse_basic(pro_euka_sequence,subline.str[subcol][z],&realnum);
            }
            pro_euka_sequence[realnum++] = '\n';
            free(subline.str);
        }

        else if (i == line_end){ //输出到col_end
            split(line," ",&subline); 
            int subcol = col_end / 10 + 1;
            int col_over = col_end % 10;
            
            int col = 1;
            for (;col < subcol;++col)
            {
                for (z = 0;z < 10;++z)
                    reverse_basic(pro_euka_sequence,subline.str[col][z],&realnum);
            }
            int k = 0;
            for (;k <= col_over;++k)
            {
                reverse_basic(pro_euka_sequence,subline.str[subcol][k],&realnum);
            }
            free(subline.str);
            break;
        }
        i++;
    }while(fgets(line,MAXCOL,fp));
    pro_euka_sequence[realnum] = '\0';
    printf("\n");
}
void reverse_basic(char group[],char basic,int *j)
{
    if (basic == 'a') group[*j] = 't';
    else if (basic == 't') group[*j] = 'a';
    else if (basic == 'c') group[*j] = 'g';
    else if (basic == 'g') group[*j] = 'c';
    else group[*j] = basic;
    (*j)++;
}
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
void Reverse(char s[])
{
     for(int i =0,j=strlen(s)-1;i<j; ++i,--j)

    {
          char c=s[i];
          s[i]=s[j];
          s[j]=c;
     }
}
#endif //_libgenbank
