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

char *data_join[SECTION] = {NULL};         //�洢join����ʼ��ֹλ����Ϣ
char *data_complement[SECTION] = {NULL};   //�洢complement����ʼ����ֹλ����Ϣ
int num_join = 0;                        //��Ҫ�������ӵ��������
int num_complement = 0;                  //���������������
char defination[MAXCOL];                 //�洢genbank��defination������
char note_fasta[MAXCOL];                 //�洢fasta�ļ�ע�������������
char complement_output[MAXCOL];
int num_array_complement = 0;

int split(char *src, char *delim, IString* istr);     //�ָ��
//void Reverse(char s[]);                               //���ú���
void tackle_join(FILE *fp,FILE *target);                //����join���͵�����
void tackle_complement(FILE *fp,FILE *target);          //����complement���͵�����
void get_sequence(int begin,int end,FILE *fp,int* i,FILE *target);       //��ȡ�����ڵ�ĸ������
void get_sequence_p(int begin,int end,FILE *fp,int *i,FILE *target);     //��ȡ�����ڵĻ���������
void reverse_basic(char basic);                                 //��������ת���ַ�
void get_CDS(char name[],FILE *target);                 //����genbank��CDS����
void get_note(FILE *fp,FILE *target);                   //����genbank��Ҫ���뵽fasta�ļ���ע��������

void get_note(FILE *fp,FILE *target)
{
    int i = 0;
    note_fasta[i++] = '>';
    int num_def = 0;
    while(fgets(line,MAXCOL,fp))
    {
        split(line," ",&subline);
        if (strcmp(subline.str[0],"DEFINITION") == 0){          //��defination�����ݿ�������
            for (int z = 1; z < subline.num; ++z){
                defination[num_def++] = ' ';
                strcpy(defination + num_def,subline.str[z]);
                num_def += strlen(subline.str[z]);
            }
            free(subline.str);
        }  
            
        else if (strcmp(subline.str[0],"VERSION") == 0){       //����version
            for (int j = 1; j < subline.num; ++j){
                if (subline.str[j][0] == 'G'){                 //����GI
                    note_fasta[i++] = 'g';
                    note_fasta[i++] = 'i';
                    note_fasta[i++] = '|';
                    strcpy(note_fasta + i,subline.str[j] + 3);
                    i += strlen(subline.str[j] + 3);
                    if (subline.str[j][strlen(subline.str[j]) - 1] == '\n') note_fasta[i - 1] = '|'; //��ֹ���һ���ַ�Ϊ\n�ƻ���ʽ
                    else note_fasta[i++] = '|';
                }
                else {                                       //�������������
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
void get_CDS(char name[],FILE *target)          //����CDSԤ����
{
    FILE *fp = fopen(name,"r");
    get_note(fp,target);
    while(fgets(line,MAXCOL,fp))               //��ȡ�ļ�ֱ��CDS��һ��
    {
        split(line," ",&subline);
        if (strcmp(subline.str[0],"CDS") == 0){
            split(subline.str[1],",<>()",&data);    //���ֳ�����ķ�ʽ�Լ��������ʼ��ֹλ�ã�����data
                                                    //����\n�Ĵ��ڣ���ˣ�data.numӦ�ü�һ
            if (data.str[0][0] >= '0' && data.str[0][0] <= '9'){
                data_join[num_join++] = data.str[0];
                break;
            }
            else
            {
                for (int i = 0; i < data.num - 1;)
                {
                    if (data.str[i][0] == 'j'){                 //����CDS��join
                        i++;
                        do
                        {
                            data_join[num_join++] = data.str[i++];
                        }while (i < data.num - 1 && data.str[i][0] >= '0' && data.str[i][0] <= '9');
                    }
                    else if (data.str[i][0] == 'c'){            //����CDS��complement
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
void tackle_complement(FILE *fp,FILE *target)      //����complement
{
    int *section_start = (int *)malloc(sizeof(int) * num_complement);
    int *section_end = (int *)malloc(sizeof(int) * num_complement);
    int j = 0;
    for (int i = 0; i < num_complement; i++)       //����ʼ��ֹ��Ƭ�δ洢��������
    {
        split(data_complement[i],".",&location);
        section_start[j] = atoi(location.str[0]);
        section_end[j++] = atoi(location.str[1]);
        free(location.str);
    }

    char sequence[] = "1";
    while(fgets(line,MAXCOL,fp))                    //��λ�������� 
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
void tackle_join(FILE *fp,FILE *target)            //����join
{
    int *section_start = (int *)malloc(sizeof(int) * num_join);
    int *section_end = (int *)malloc(sizeof(int) * num_join);
    int j = 0;
    for (int i = 0; i < num_join; i++)            //����ʼ��ֹ��Ƭ�δ洢��������
    {
        split(data_join[i],".",&location);
        section_start[j] = atoi(location.str[0]);
        section_end[j++] = atoi(location.str[1]);
        free(location.str);
    }
    
    char sequence[] = "1";
    while(fgets(line,MAXCOL,fp))  //���ļ����ת�Ƶ�������ʼ����һ��
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
    fseek(fp,0,SEEK_SET); //����join��complement�������������Ҫ�ض���
    free(section_start);
    free(section_end);
}

void get_sequence(int begin,int end,FILE *fp,int* i,FILE *target)    //����join���ж�ȡ
{
    int line_begin = (begin - 1) / 60;
    int line_end = (end - 1) / 60;
    int col_begin = (begin - 1) % 60;
    int col_end = (end - 1) % 60;
    int z;
    do  //���ж�ȡ
    {
        if ((*i) == line_begin){
            split(line," ",&subline); //�������
            int subcol = col_begin / 10 + 1; //��ȡ�ִ����ڵ�subline��λ��
            int col_start = col_begin % 10;  //��ȡ����
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

        else if ((*i) > line_begin && (*i) < line_end){ //�������
            int subcol = 1;
            split(line," ",&subline); //�������
            for (;subcol < 7;++subcol)
            {
                for (z = 0;z < 10;++z)
                    writeFasta_seq(subline.str[subcol][z],target);
            }
            free(subline.str);
        }

        else if ((*i) == line_end){ //�����col_end
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
void get_sequence_p(int begin,int end,FILE *fp,int *i,FILE *target) //����complement���ж�ȡ
{
    int line_begin = (begin - 1) / 60;
    int line_end = (end - 1) / 60;
    int col_begin = (begin - 1) % 60;
    int col_end = (end - 1) % 60;

    int z; //��������
    do  //���ж�ȡ
    {
        if ((*i) == line_begin){
            num_array_complement = 0;
            split(line," ",&subline); //�������
            int subcol = col_begin / 10 + 1; //��ȡ�ִ����ڵ�subline��λ��
            int col_start = col_begin % 10;  //��ȡ����
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
            //complement_sequence[realnum++] = '\n';    //ʹ����������ֻ�����в���
            free(subline.str);
        }

        else if ((*i) > line_begin && (*i) < line_end){ //�������
            int subcol = 1;
            split(line," ",&subline); //�������
            for (;subcol < 7;++subcol)
            {
                for (z = 0;z < 10;++z)
                    reverse_basic(subline.str[subcol][z]);
            }
            free(subline.str);
        }

        else if ((*i) == line_end){ //�����col_end
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
int split(char *src, char *delim, IString* istr)     //split src�����ָ����Ӵ�����ṹ���str��
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