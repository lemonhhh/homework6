/*
 * @Descripttion: 
 * @version: 
 * @Author: sueRimn
 * @Date: 2020-06-11 20:23:47
 * @LastEditors: sueRimn
 * @LastEditTime: 2020-06-11 20:24:19
 */ 
#include"libfasta.h"

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