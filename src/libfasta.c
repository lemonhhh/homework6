#include "include/libfasta.h"

void writeFasta_seq(int c, FILE *fp)
{
    if (currentcol < 70)
    { //������г��� + ԭ�е�С��70��ֱ������
        putc(c, fp);
        currentcol += 1;
    }
    else
    {
        putc('\n', fp); //�����һ��
        putc(c, fp);
        currentcol = 1;
    }
}

void writeFasta_note(char note[], FILE *fp) //�Դ����з�������ע����Ϣ
{
    fputs(note, fp);
}