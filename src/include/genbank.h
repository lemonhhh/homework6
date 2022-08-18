/*
 * @Descripttion: 
 * @version: 
 * @Author: sueRimn
 * @Date: 2020-06-05 11:39:49
 * @LastEditors: sueRimn
 * @LastEditTime: 2020-06-11 20:24:27
 */ 
#ifndef _libgenbank
#define _libgenbank

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

#define MAXCOL 5000
#define SECTION 100

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
int *section_start_c;
int *section_end_c;
int *section_start_j;
int *section_end_j;

int split(char *src, char *delim, IString* istr);     //分割函数
//void Reverse(char s[]);                               //倒置函数
void tackle_join(FILE *fp,FILE *target);                //处理join类型的连接
void tackle_complement(FILE *fp,FILE *target);          //处理complement类型的连接
void get_sequence(int begin,int end,FILE *fp,int* i,FILE *target);       //获取区间内的母链序列
void get_sequence_p(int begin,int end,FILE *fp,int *i,FILE *target);     //获取区间内的互补链序列
void reverse_basic(char basic);                                 //互补链，转置字符
void get_CDS(char name[],FILE *target);                 //处理genbank的CDS内容
void get_note(FILE *fp,FILE *target);                   //处理genbank需要输入到fasta文件的注释行内容
void output(FILE *fp,FILE *target);




#endif //_libgenbank