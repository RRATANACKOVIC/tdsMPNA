#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void saxpy(int size, float a, float* x, float* y)
{
	#pragma omp parallel for
	for(int i = 0; i<size; i++)
	{
		y[i] = a * y[i] + x[i];
	}
}

float dotprod (int size, float beta, float* x, float* y)
{
	float output = beta;
	#pragma omp parallel reduction(+:output)	
	for(int i = 0; i<size; i++)
	{
		output += y[i] * x[i];
	}
	return output;
}

float norm (int size, float* x)
{
		return	sqrt(dotprod(size, 0.0, x, x));
}

void printvec(int size, float* input)
{
	printf("vector display (size = %d)\n", size);
	
	for(int i = 0; i<size; i++)
	{
		printf("vector[%d] = %f\n",i,*(input+i));
	}
}
int main(void)
{
	float x[5] = { 1.0, 2.0, 3.0, 4.0, 5.0};// 15, rac(15) ~ 4
	float y[5] = { 6.0, 7.0, 8.0, 9.0, 10.0};// 40, rac (40) ~ 6
	//              6    14   24   36   50 = 20 + 60 + 50 = 130 
	float beta = 10.0;
	printf("dotprod %d + x^T.y = %f\n", (int)beta, dotprod(5, beta,x,y)); 
	printf("||x|| = %f\n", norm(5,x)); 
	printf("||y|| = %f\n", norm(5,y)); 
	saxpy(5, 1.0, x, y);
	printvec(5, y);
	return 0;
}
