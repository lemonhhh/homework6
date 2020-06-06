/*
 * @Descripttion: 
 * @version: 
 * @Author: sueRimn
 * @Date: 2020-06-05 10:16:32
 * @LastEditors: sueRimn
 * @LastEditTime: 2020-06-06 11:42:20
 */ 
#ifndef _fasta_h
#define _fasta_h

#include<stdio.h>
#include<stdlib.h>
#include<string.h>

int currentcol = 0;

void writeFasta_seq(int c,FILE *fp);        //写入fasta文件的序列部分   
void writeFasta_note(char note[],FILE *fp);         

void writeFasta_seq(int c,FILE *fp)
{
    if (currentcol < 70){   //如果序列长度 + 原有的小于70，直接输入
        putc(c,fp);
        currentcol += 1;
    }
    else{
        putc('\n',fp);   //光标下一行
        putc(c,fp);
        currentcol = 1;
    }
}
void writeFasta_note(char note[],FILE *fp)   //自带换行符，输入注释信息
{
    fputs(note,fp);
}
#endif  //_fasta_h