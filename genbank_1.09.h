/*
 * @Descripttion: 
 * @version: 
 * @Author: sueRimn
 * @Date: 2020-06-05 11:39:49
 * @LastEditors: sueRimn
 * @LastEditTime: 2020-06-08 20:50:56
 */ 
#ifndef _libgenbank
#define _libgenbank

#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include "libfasta.h"

#define MAXCOL 2000
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

char *data_join[SECTION] = {NULL};         //存储join的起始终止位置信息
char *data_complement[SECTION] = {NULL};   //存储complement的起始和终止位置信息
int num_join = 0;                        //需要正向连接的区间个数
int num_complement = 0;                  //互补链的区间个数
char defination[MAXCOL];                 //存储genbank的defination的内容
char note_fasta[MAXCOL];                 //存储fasta文件注释行输入的内容
char complement_output[MAXCOL];
int num_array_complement = 0;

int split(char *src, char *delim, IString* istr);     //分割函数
//void Reverse(char s[]);                               //倒置函数
void tackle_join(FILE *fp,FILE *target);                //处理join类型的连接
void tackle_complement(FILE *fp,FILE *target);          //处理complement类型的连接
void get_sequence(int begin,int end,FILE *fp,int* i,FILE *target);       //获取区间内的母链序列
void get_sequence_p(int begin,int end,FILE *fp,int *i,FILE *target);     //获取区间内的互补链序列
void reverse_basic(char basic);                                 //互补链，转置字符
void get_CDS(char name[],FILE *target);                 //处理genbank的CDS内容
void get_note(FILE *fp,FILE *target);                   //处理genbank需要输入到fasta文件的注释行内容

void get_note(FILE *fp,FILE *target)
{
    int i = 0;
    note_fasta[i++] = '>';
    int num_def = 0;
    while(fgets(line,MAXCOL,fp))
    {
        split(line," ",&subline);
        if (strcmp(subline.str[0],"DEFINITION") == 0){          //把defination的数据拷贝下来
            for (int z = 1; z < subline.num; ++z){
                defination[num_def++] = ' ';
                strcpy(defination + num_def,subline.str[z]);
                num_def += strlen(subline.str[z]);
            }
            free(subline.str);
        }  
            
        else if (strcmp(subline.str[0],"VERSION") == 0){       //处理version
            for (int j = 1; j < subline.num; ++j){
                if (subline.str[j][0] == 'G'){                 //处理GI
                    note_fasta[i++] = 'g';
                    note_fasta[i++] = 'i';
                    note_fasta[i++] = '|';
                    strcpy(note_fasta + i,subline.str[j] + 3);
                    i += strlen(subline.str[j] + 3);
                    if (subline.str[j][strlen(subline.str[j]) - 1] == '\n') note_fasta[i - 1] = '|'; //防止最后一个字符为\n破坏格式
                    else note_fasta[i++] = '|';
                }
                else {                                       //其他的引用情况
                    note_fasta[i++] = 'r';
                    note_fasta[i++] = 'e';
                    note_fasta[i++] = 'f';
                    note_fasta[i++] = '|';
                    strcpy(note_fasta + i,subline.str[j]);
                    i += strlen(subline.str[j]);
                    if (subline.str[j][strlen(subline.str[j]) - 1] == '\n') note_fasta[i - 1] = '|';
                    else note_fasta[i++] = '|';
                }
            }
        }
        else if (strcmp(subline.str[0],"DBLINK") == 0) 
            break;
    }
    strcpy(note_fasta + i,defination);
    writeFasta_note(note_fasta,target);
}
void get_CDS(char name[],FILE *target)          //进行CDS预处理
{
    FILE *fp = fopen(name,"r");
    get_note(fp,target);
    while(fgets(line,MAXCOL,fp))               //读取文件直至CDS那一行
    {
        split(line," ",&subline);
        if (strcmp(subline.str[0],"CDS") == 0){
            split(subline.str[1],",<>()",&data);    //划分出编码的方式以及编码的起始终止位置，存在data
                                                    //由于\n的存在，因此，data.num应该减一
            if (data.str[0][0] >= '0' && data.str[0][0] <= '9'){
                data_join[num_join++] = data.str[0];
                break;
            }
            else
            {
                for (int i = 0; i < data.num - 1;)
                {
                    if (data.str[i][0] == 'j'){                 //处理CDS的join
                        i++;
                        do
                        {
                            data_join[num_join++] = data.str[i++];
                        }while (i < data.num - 1 && data.str[i][0] >= '0' && data.str[i][0] <= '9');
                    }
                    else if (data.str[i][0] == 'c'){            //处理CDS的complement
                        i++;
                        do
                        {
                            data_complement[num_complement++] = data.str[i++];
                        }while (i < data.num - 1 && data.str[i][0] >= '0' && data.str[i][0] <= '9' );
                    }
                }
            }
            
        }
        else free(subline.str);
    }
    if (num_join != 0)
        tackle_join(fp,target);
    if (num_complement != 0)
        tackle_complement(fp,target);
    char end[] = "//";
    writeFasta_note(end,target);
    fclose(fp);
}
void tackle_complement(FILE *fp,FILE *target)      //处理complement
{
    int *section_start = (int *)malloc(sizeof(int) * num_complement);
    int *section_end = (int *)malloc(sizeof(int) * num_complement);
    int j = 0;
    for (int i = 0; i < num_complement; i++)       //把起始终止的片段存储在数组中
    {
        split(data_complement[i],".",&location);
        section_start[j] = atoi(location.str[0]);
        section_end[j++] = atoi(location.str[1]);
        free(location.str);
    }

    char sequence[] = "1";
    while(fgets(line,MAXCOL,fp))                    //定位到序列行 
    {
        split(line," ",&subline);
        if (strcmp(subline.str[0],sequence) == 0)
            break;
        else free(subline.str);
    }
    int line_num_ = 0;
    for (int k = 0; k < num_complement;++k)
    {
        //printf("Section %d:\n",k + 1);
        get_sequence_p(section_start[k],section_end[k],fp,&line_num_,target);
        //printf("%s\n",complement_sequence);
    }
    //Reverse(complement_sequence);
    writeFasta_seq('\n',target);
    free(section_start);
    free(section_end);
}
void tackle_join(FILE *fp,FILE *target)            //处理join
{
    int *section_start = (int *)malloc(sizeof(int) * num_join);
    int *section_end = (int *)malloc(sizeof(int) * num_join);
    int j = 0;
    for (int i = 0; i < num_join; i++)            //把起始终止的片段存储在数组中
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
    for (int k = 0; k < num_join;++k)
    {
        //printf("Section %d:\n",k + 1);
        get_sequence(section_start[k],section_end[k],fp,&line_num_,target);
        //writeFasta(complement_sequence,target);
    }
    fseek(fp,0,SEEK_SET); //由于join与complement并不连续，因此要重定向
    free(section_start);
    free(section_end);
}

void get_sequence(int begin,int end,FILE *fp,int* i,FILE *target)    //处理join序列读取
{
    int line_begin = (begin - 1) / 60;
    int line_end = (end - 1) / 60;
    int col_begin = (begin - 1) % 60;
    int col_end = (end - 1) % 60;
    int z;
    do  //序列读取
    {
        if ((*i) == line_begin){
            split(line," ",&subline); //读到最后
            int subcol = col_begin / 10 + 1; //获取字串所在的subline的位置
            int col_start = col_begin % 10;  //获取列数
            for (;col_start < 10;++col_start)
            {
                writeFasta_seq(subline.str[subcol][col_start],target);
            }
            subcol++;
            for (;subcol < 7;++subcol)
            {
                for (z = 0;z < 10;++z)
                    writeFasta_seq(subline.str[subcol][z],target);
            }
            free(subline.str);
        }

        else if ((*i) > line_begin && (*i) < line_end){ //输出整行
            int subcol = 1;
            split(line," ",&subline); //读到最后
            for (;subcol < 7;++subcol)
            {
                for (z = 0;z < 10;++z)
                    writeFasta_seq(subline.str[subcol][z],target);
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
                for (z = 0;z < 10;++z)
                    writeFasta_seq(subline.str[col][z],target);
            }
            int k = 0;
            for (;k <= col_over;++k)
            {
                writeFasta_seq(subline.str[subcol][k],target);
            }
            free(subline.str);
            break;
        }
        
        (*i)++;
    }while(fgets(line,MAXCOL,fp));

}
void get_sequence_p(int begin,int end,FILE *fp,int *i,FILE *target) //处理complement序列读取
{
    int line_begin = (begin - 1) / 60;
    int line_end = (end - 1) / 60;
    int col_begin = (begin - 1) % 60;
    int col_end = (end - 1) % 60;

    int z; //计数因子
    do  //序列读取
    {
        if ((*i) == line_begin){
            num_array_complement = 0;
            split(line," ",&subline); //读到最后
            int subcol = col_begin / 10 + 1; //获取字串所在的subline的位置
            int col_start = col_begin % 10;  //获取列数
            for (;col_start < 10;++col_start)
            {
                reverse_basic(subline.str[subcol][col_start]);
            }
            subcol++;
            for (;subcol < 7;++subcol)
            {
                for (z = 0;z < 10;++z)
                    reverse_basic(subline.str[subcol][z]);
            }
            //complement_sequence[realnum++] = '\n';    //使整个数组中只有序列部分
            free(subline.str);
        }

        else if ((*i) > line_begin && (*i) < line_end){ //输出整行
            int subcol = 1;
            split(line," ",&subline); //读到最后
            for (;subcol < 7;++subcol)
            {
                for (z = 0;z < 10;++z)
                    reverse_basic(subline.str[subcol][z]);
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
                for (z = 0;z < 10;++z)
                    reverse_basic(subline.str[col][z]);
            }
            int k = 0;
            for (;k <= col_over;++k)
            {
                reverse_basic(subline.str[subcol][k]);
            }
            free(subline.str);
            for (int i = num_array_complement - 1; i >= 0; --i)
                writeFasta_seq(complement_output[i],target);
            break;
        }
        (*i)++;
    }while(fgets(line,MAXCOL,fp));
}
void reverse_basic(char basic)
{
    if (basic == 'a') complement_output[num_array_complement++] = 't';
    else if (basic == 't') complement_output[num_array_complement++] = 'a';
    else if (basic == 'c') complement_output[num_array_complement++] = 'g';
    else if (basic == 'g') complement_output[num_array_complement++] = 'c';
}
int split(char *src, char *delim, IString* istr)     //split src，将分割后的子串存入结构体的str中
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
/*void Reverse(char s[])
{
     for(int i =0,j=strlen(s)-1;i<j; ++i,--j)

    {
          char c=s[i];
          s[i]=s[j];
          s[j]=c;
     }
}*/


#endif //_libgenbank