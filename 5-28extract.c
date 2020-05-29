#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

#define MAXLEN 1000
#define MAX 10
#define MAXNUM 10000

int main(int argc,char **argv)
{
   FILE *fp;
   FILE *fa;

   int flag=0;
   int start,end;
   char Start[MAX];
   char End[MAX];
   char line[MAXLEN];
   char gene[MAXNUM];

   int i=0,j=0,k=0;
   int count=0;

   fp=fopen(argv[1],"r");

   while(fgets(line,MAXLEN,fp))
   {
     if(strstr(line,"CDS"))
     {

        if(strstr(line,"complement"))
        {
          flag=1;

        }

        while(!isdigit(line[i])) i++;
        while(isdigit(line[i]))
        {
            Start[j++]=line[i++];
        }
        Start[j]='\0';
        //printf("%s\n",Start);
        while(!isdigit(line[i])) i++;
        while(isdigit(line[i]))
        {
            End[k++]=line[i++];
        }
        End[k]='\0';
        //printf("%s\n",End);
        start=atoi(Start);
        end=atoi(End);
        printf("length is %d\n",end-start+1);


     }

     if(strstr(line,"ORIGIN"))
     {
       //printf("%s\n",line);
       break;
     }

   }


   j=0;

   while(fgets(line,MAXLEN,fp))
   {
     i=0;
      while(!isalpha(line[i])) i++;
      //printf("%c\n",line[i]);

       while(line[i])
      {
        if(isalpha(line[i]))
        {
          if(count<start)
          {
            count++;
          }
          if(count>=start && count<=end)
          {
            gene[j]=line[i];
            count++;
            j++;
          }


        }
         i++;
      }


   }

   gene[j]='\0';

   
   printf("%s\n",gene);
   // printf("%ld\n",strlen(gene));


   fa=fopen("output.fasta","w+");
   fputs(gene,fa);
   fclose(fa);

   return 0;
}
