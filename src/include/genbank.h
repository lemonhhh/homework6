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

char *data_join[SECTION] = {NULL};         //�洢join����ʼ��ֹλ����Ϣ
char *data_complement[SECTION] = {NULL};   //�洢complement����ʼ����ֹλ����Ϣ
int num_join = 0;                        //��Ҫ�������ӵ��������
int num_complement = 0;                  //���������������
char defination[MAXCOL];                 //�洢genbank��defination������
char note_fasta[MAXCOL];                 //�洢fasta�ļ�ע�������������
char complement_output[MAXCOL];
int num_array_complement = 0;
int *section_start_c;
int *section_end_c;
int *section_start_j;
int *section_end_j;

int split(char *src, char *delim, IString* istr);     //�ָ��
//void Reverse(char s[]);                               //���ú���
void tackle_join(FILE *fp,FILE *target);                //����join���͵�����
void tackle_complement(FILE *fp,FILE *target);          //����complement���͵�����
void get_sequence(int begin,int end,FILE *fp,int* i,FILE *target);       //��ȡ�����ڵ�ĸ������
void get_sequence_p(int begin,int end,FILE *fp,int *i,FILE *target);     //��ȡ�����ڵĻ���������
void reverse_basic(char basic);                                 //��������ת���ַ�
void get_CDS(char name[],FILE *target);                 //����genbank��CDS����
void get_note(FILE *fp,FILE *target);                   //����genbank��Ҫ���뵽fasta�ļ���ע��������
void output(FILE *fp,FILE *target);




#endif //_libgenbank