/***
	sort two array by the order of the first one
	Qsort is implemented following the algorithm explained in Numrec C, pp 332-335
	Under a certain size, a part of the array is sorted by the naive method (straight insertion)
***/

#include <stdlib.h>
#include <stdio.h>

#define NSTACK 1000                           /* the number of subarray that could be memorized */
#define SWITCH_METHOD 7                       /* under that size switch to a Insertion method */

#define SWAPf(a,b) {tmpf=(a); (a)=(b); (b)=tmpf;}
#define SWAPi(a,b) {tmpi=(a); (a)=(b); (b)=tmpi;}

void qsort2(float *array1, long *array2, long n){

	
	float *Stack;       /* where remaining values of subarrays are temporarily stocked */
	long nStack=0;     /* How many of these values are in the Stack. Is never odd */
	
	long end=n-1,      /* the last value of the array to sort */
	     beg=0,        /* the first value ---  */
		  postbeg;      /* the second value ---  */

	float val_postbeg;  /* the value of the postbeg position - the one which is used for partion-exchange */

	long demi_len;     /* half of the distance between beg and end */
	
	long i,            /* counter from postbeg to right */
	     j;            /* counter from end to left */

	float val1_i,       /* for insertion stock temporarily a value */ 
              val2_i;
	
	float tmpf;          /* used for the SWAPf macro */
	long tmpi;          /* used for the SWAPi macro */
	
	
	Stack = (float *)malloc( (size_t) NSTACK*sizeof(float));
	if(! Stack )fprintf(stderr , "qsort2: not enough memory for Stack"), exit(1) ; 
		
	while(1){ 

		if( end-beg+1 > SWITCH_METHOD){
	
			demi_len = (end-beg) >> 1 ; 
			postbeg = beg+1;
			
			SWAPf( array1[beg+demi_len], array1[postbeg] );
			SWAPi( array2[beg+demi_len], array2[postbeg] );
			
			if(array1[beg] > array1[postbeg]){            /* rearrange to have  beg <= postbeg <= end */ 
				SWAPf( array1[beg], array1[postbeg] );
				SWAPi( array2[beg], array2[postbeg] );
			}

			if(array1[beg] > array1[end]){
				SWAPf( array1[beg], array1[end] );
				SWAPi( array2[beg], array2[end] );
			}

			if(array1[postbeg] > array1[end]){
				SWAPf( array1[postbeg], array1[end] );
				SWAPi( array2[postbeg], array2[end] );
			}
			
			
			i = postbeg;
			j = end;
						
			val_postbeg =  array1[postbeg];
			
			while(1)                                   /* enter the partition exchange process */
				{
					do i++; while( array1[i] < val_postbeg );
					do j--; while( array1[j] > val_postbeg );
					
					if(j<i) break;
					
					SWAPf( array1[i], array1[j] );
					SWAPi( array2[i], array2[j] );
				}
						
			SWAPf( array1[postbeg] , array1[j] );   /* place the postbeg value into j */
			SWAPi( array2[postbeg] , array2[j] );
			
			if(nStack+2 > NSTACK)
				fprintf(stderr, "qsort2: not enough Stack... sorry bye\n"),exit(1);
			
			if(end-i+1 >= j-beg){
				Stack[nStack++] = i;                /* stock in Stack the largest and go with the smallest */
				Stack[nStack++] = end;
				
				end = j-1;
			}
			else{
				Stack[nStack++] = beg;
				Stack[nStack++] = j-1;
				
				beg = i;
			}
			
		}
		else{                                      /* Under a certain size, switch to the straight insertion method */
	
			
			for(i=beg+1; i <= end ; i++){
			
				val1_i = array1[i];
				val2_i = array2[i];
				
				for(j=i-1;j>=beg;j--){
				
					if(array1[j] < val1_i)break;
					
					array1[j+1] = array1[j];
					array2[j+1] = array2[j];
				}
				array1[j+1]=val1_i;
				array2[j+1]=val2_i;
			}

			if(nStack==0)break;          /* this si the end  - exit the process */
			
			end = Stack[--nStack];       /* go for the next round with the stacked parameters */
			beg = Stack[--nStack];
		}

	 }

	free(Stack);

}
