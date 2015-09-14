#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<iostream>
#include <ctype.h>
#include"misc.h"

using namespace std;

/*The following function check if the input string is a valid RNA sequence*/
int CheckSequence(char * sequence){     //safe
  int i;
  for (i=1;i<strlen(sequence);++i)
    {
     if (toupper(sequence[i])!='A' && toupper(sequence[i])!='U' && toupper(sequence[i])!='C' && toupper(sequence[i])!='G'&& toupper(sequence[i])!='T') //check if there are invalid characters
        {
         printf("The input string should only contain A,U,T,C,G!\n");
         exit(1);
        }
     else if (toupper(sequence[i])=='T') //change T to U
        sequence[i]='U';
     else
        sequence[i]=toupper(sequence[i]); // change lower case to upper case
    }
  return 0;
}
