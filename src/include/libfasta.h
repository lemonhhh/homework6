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


#endif  //_fasta_h