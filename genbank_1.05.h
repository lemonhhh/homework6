#ifndef _libgenbank
#define _libgenbank

#include<stdlib.h>
#include<string.h>
#include<stdio.h>

#define MAXCOL 1000
#define SECTION 10
typedef struct {
    char **str;     //the PChar of string array
    size_t num;     //the number of string
}IString;

char buf[MAXCOL];
char line[MAXCOL];
IString subline;
IString data;
IString location;
char *complement_sequence; //存储互补链的碱基
char *data_join[SECTION] = {NULL};         //存储join的起始终止位置信息
char *data_complement[SECTION] = {NULL};   //存储complement的起始和终止位置信息
int num_join = 0;
int num_complement = 0;

int split(char *src, char *delim, IString* istr);     //分割函数
void Reverse(char s[]);                               //倒置函数
void tackle_join(FILE *fp);
void tackle_complement(FILE *fp);
void get_sequence(int begin,int end,FILE *fp,int* i);
void get_sequence_p(int begin,int end,FILE *fp,int *i);
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
int main()
{
    get_CDS("c:/vsproject/zhan_lian/Linux.txt");
    system("PAUSE");
    return 0;
}
void get_CDS(char name[])          //进行CDS预处理
{
    FILE *fp = fopen(name,"r");
    while(fgets(line,MAXCOL,fp))   //读取文件直至CDS那一行
    {
        split(line," ",&subline);
        if (strcmp(subline.str[0],"CDS") == 0)
            break;
        else free(subline.str);
    }
    split(subline.str[1],",<>()",&data);    //划分出编码的方式以及编码的起始终止位置，存在data
                                       //由于\n的存在，因此，data.num应该减一
    for (int i = 0; i < data.num - 1;)
    {
        if (data.str[i][0] == 'j'){
            i++;
            do
            {
                data_join[num_join++] = data.str[i++];
            }while (i < data.num - 1 && data.str[i][0] >= '0' && data.str[i][0] <= '9');
        }
        else if (data.str[i][0] == 'c'){
            i++;
            do
            {
                data_complement[num_complement++] = data.str[i++];
            }while (i < data.num - 1 && data.str[i][0] >= '0' && data.str[i][0] <= '9' );
        }
    }
    free(data.str);
    if (num_join != 0)
        tackle_join(fp);
    if (num_complement != 0)
        tackle_complement(fp);

    fclose(fp);
}
void tackle_complement(FILE *fp)    //处理complement
{
    int *section_start = (int *)malloc(sizeof(int) * num_complement);
    int *section_end = (int *)malloc(sizeof(int) * num_complement);
    int j = 0;
    for (int i = 0; i < num_complement; i++)  //把起始终止的片段存储在数组中
    {
        split(data_complement[i],".",&location);
        section_start[j] = atoi(location.str[0]);
        section_end[j++] = atoi(location.str[1]);
        free(location.str);
    }

    char sequence[] = "1";
    while(fgets(line,MAXCOL,fp)) //定位到序列行 
    {
        split(line," ",&subline);
        if (strcmp(subline.str[0],sequence) == 0)
            break;
        else free(subline.str);
    }
    int line_num_ = 0;
    printf("位于互补链的序列为：\n");
    for (int k = 0; k < num_complement;++k)
    {
        printf("Section %d:\n",k + 1);
        get_sequence_p(section_start[k],section_end[k],fp,&line_num_);
        printf("%s\n",complement_sequence);
    }
    //Reverse(complement_sequence);
    free(section_start);
    free(section_end);
}
void tackle_join(FILE *fp)      //处理join
{
    int *section_start = (int *)malloc(sizeof(int) * num_join);
    int *section_end = (int *)malloc(sizeof(int) * num_join);
    int j = 0;
    for (int i = 0; i < num_join; i++)  //把起始终止的片段存储在数组中
    {
        split(data_join[i],".",&location);
        section_start[j] = atoi(location.str[0]);
        section_end[j++] = atoi(location.str[1]);
        free(location.str);
    }
    
    char sequence[] = "1";
    while(fgets(line,MAXCOL,fp))  //将文件光标转移到序列起始的那一行
    {
        split(line," ",&subline);
        if (strcmp(subline.str[0],sequence) == 0)
            break;
        else free(subline.str);
    }
    int line_num_ = 0;
    printf("位于母链的序列为：\n");
    for (int k = 0; k < num_join;++k)
    {
        printf("Section %d:\n",k + 1);
        get_sequence(section_start[k],section_end[k],fp,&line_num_);
    }
    fseek(fp,0,SEEK_SET); //由于join与complement并不连续，因此要重定向
    free(section_start);
    free(section_end);
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
void get_sequence_p(int begin,int end,FILE *fp,int *i) //原核生物基因 CDS序列输出
{
    int line_begin = (begin - 1) / 60;
    int line_end = (end - 1) / 60;
    int col_begin = (begin - 1) % 60;
    int col_end = (end - 1) % 60;
    free(complement_sequence);
    complement_sequence = (char *)malloc(sizeof(char) * (end - begin + 300)); //将序列存储到字符串中(还有换行符)
    int realnum = 0;

    int z; //计数因子
    do  //序列读取
    {
        if ((*i) == line_begin){
            split(line," ",&subline); //读到最后
            int subcol = col_begin / 10 + 1; //获取字串所在的subline的位置
            int col_start = col_begin % 10;  //获取列数
            for (;col_start < 10;++col_start)
            {
                reverse_basic(complement_sequence,subline.str[subcol][col_start],&realnum);
            }
            subcol++;
            for (;subcol < 7;++subcol)
            {
                for (z = 0;z < 10;++z)
                    reverse_basic(complement_sequence,subline.str[subcol][z],&realnum);
            }
            complement_sequence[realnum++] = '\n';
            free(subline.str);
        }

        else if ((*i) > line_begin && (*i) < line_end){ //输出整行
            int subcol = 1;
            split(line," ",&subline); //读到最后
            for (;subcol < 7;++subcol)
            {
                for (z = 0;z < 10;++z)
                    reverse_basic(complement_sequence,subline.str[subcol][z],&realnum);
            }
            complement_sequence[realnum++] = '\n';
            free(subline.str);
        }

        else if ((*i) == line_end){ //输出到col_end
            split(line," ",&subline); 
            int subcol = col_end / 10 + 1;
            int col_over = col_end % 10;
            
            int col = 1;
            for (;col < subcol;++col)
            {
                for (z = 0;z < 10;++z)
                    reverse_basic(complement_sequence,subline.str[col][z],&realnum);
            }
            int k = 0;
            for (;k <= col_over;++k)
            {
                reverse_basic(complement_sequence,subline.str[subcol][k],&realnum);
            }
            free(subline.str);
            break;
        }
        (*i)++;
    }while(fgets(line,MAXCOL,fp));
    complement_sequence[realnum] = '\0';
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
