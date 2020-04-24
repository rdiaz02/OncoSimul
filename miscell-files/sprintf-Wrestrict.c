#include <stdio.h>
#include <string.h>
#include <math.h>

int main () {
   char str[8000];

   sprintf(str, "Value of Pi = %f", M_PI);
   //sprintf(str, "%s y coco = %f", str, 95.5);
   sprintf(str + strlen(str), "y coco = %f",  95.5);
   
   //sprintf(str, "%s otras cosas", str);

   
   puts(str);
   
   return(0);
}
