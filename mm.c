#include<stdio.h>
#include<stdlib.h>
#include<time.h>
#include "mm.h"



// Task 1: Flush the cache so that we can do our measurement :)
void flush_all_caches()
{
	// Your code goes here
}

void load_matrix_base()
{
	long i;
	huge_matrixA = malloc(sizeof(long)*(long)SIZEX*(long)SIZEY);
	huge_matrixB = malloc(sizeof(long)*(long)SIZEX*(long)SIZEY);
	huge_matrixC = malloc(sizeof(long)*(long)SIZEX*(long)SIZEY);
	// Load the input
	// Note: This is suboptimal because each of these loads can be done in parallel.
	for(i=0;i<((long)SIZEX*(long)SIZEY);i++)
	{
		fscanf(fin1,"%ld", (huge_matrixA+i)); 		
		fscanf(fin2,"%ld", (huge_matrixB+i)); 		
		huge_matrixC[i] = 0;		
	}
}

void free_all()
{
	free(huge_matrixA);
	free(huge_matrixB);
	free(huge_matrixC);
}

void multiply_base()
{
	// Your code here
	//
	// Implement your baseline matrix multiply here.
	int i, j, k; 
	for(i = 0; i < SIZEX; i++){
		for(j=0 ; j < SIZEY; j++){
			for(k=0; k < SIZEX; k++){
				//multiplication on each index 
				//whereas [0*22000+0= 0    0*22000+1=1         ]
				//        [1*22000+0=22000 1*22000+1= 22001    ]
				//something along that line
				huge_matrixC[i*SIZEY+j] +=(
					huge_matrixA[i*SIZEY+k] * 
					huge_matrixB[k*SIZEY+j]
					);
			}
		}
	}

}

void compare_results()
{
	fout = fopen("./out.in","r");
	long i;
	long temp1, temp2;
	for(i=0;i<((long)SIZEX*(long)SIZEY);i++)
	{
		fscanf(fout, "%ld", &temp1);
		fscanf(ftest, "%ld", &temp2);
		if(temp1!=temp2)
		{
			printf("Wrong solution!\n");
			exit(1);
		}
	}
	printf("Right solution\n");
	fclose(fout);
	fclose(ftest);
}

void write_results()
{
	// Your code here
	//
	// Basically, make sure the result is written on fout
	// Each line represent value in the X-dimension of your matrix
	int i, j;
	for(i = 0; i < SIZEX; i++){
		for(j=0 ; j < SIZEY; j++){
			fprintf(fout, "%ld ", huge_matrixC[i*SIZEY+j]);
		}
		fprintf(fout, "\n");
	}
}

void load_matrix()
{
	// Your code here
}



void multiply()
{
	// Your code here
    //outer loop  withblocks 
	int i,j; //these are for outer loop
	int m,n; // inner loop
	int k;

	int sum_tl, sum_tr, sum_bl, sum_br;

	for(i = 0; i < SIZEX; i+=BLOCKSIZE){ //x axis of block iteration
		for(j = 0; j < SIZEY; j+=BLOCKSIZE){// y axis of block iteration
			for(m = i; m < i+BLOCKSIZE; m+=2){ 
				for(n = j; n < j+BLOCKSIZE; n+=2){
					sum_tl = 0;
					sum_tr = 0;
					sum_bl = 0;
					sum_br = 0;
					for(k=0; k < SIZEY; k++){
						sum_tl += huge_matrixA[m*SIZEY+k] * huge_matrixB[k*SIZEY+n];
						sum_tr += huge_matrixA[m*SIZEY+k] * huge_matrixB[k*SIZEY+n+1];
						sum_bl += huge_matrixA[(m+1)*SIZEY+k] * huge_matrixB[k*SIZEY+n];
						sum_br += huge_matrixA[(m+1)*SIZEY+k] * huge_matrixB[k*SIZEY+(n+1)];
					}
					huge_matrixC[m*SIZEY+n] = sum_tl;
					huge_matrixC[m*SIZEY+n+1] = sum_tr;
					huge_matrixC[(m+1)*SIZEY+n] = sum_bl;
					huge_matrixC[(m+1)*SIZEY+(n+1)] = sum_br;
				}
			}
		}
	}
}

int main()
{
	
	clock_t s,t;
	double total_in_base = 0.0;
	double total_in_your = 0.0;
	double total_mul_base = 0.0;
	double total_mul_your = 0.0;
	fin1 = fopen("./input1.in","r");
	fin2 = fopen("./input2.in","r");
	fout = fopen("./out.in","w");
	ftest = fopen("./reference.in","r");
	

	flush_all_caches();

	s = clock();
	load_matrix_base();
	t = clock();
	total_in_base += ((double)t-(double)s) / CLOCKS_PER_SEC;
	printf("[Baseline] Total time taken during the load = %f seconds\n", total_in_base);

	// s = clock();
	// multiply_base();
	// t = clock();
	// total_mul_base += ((double)t-(double)s) / CLOCKS_PER_SEC;
	// printf("[Baseline] Total time taken during the multiply = %f seconds\n", total_mul_base);
	// fclose(fin1);
	// fclose(fin2);
	// fclose(fout);
	// free_all();

	// flush_all_caches();

	// s = clock();
	// load_matrix();
	// t = clock();
	// total_in_your += ((double)t-(double)s) / CLOCKS_PER_SEC;
	// printf("Total time taken during the load = %f seconds\n", total_in_your);

	s = clock();
	multiply();
	t = clock();
	total_mul_your += ((double)t-(double)s) / CLOCKS_PER_SEC;
	printf("Total time taken during the multiply = %f seconds\n", total_mul_your);
	write_results();
	fclose(fin1);
	fclose(fin2);
	fclose(fout);
	free_all();
	compare_results();

	return 0;
}