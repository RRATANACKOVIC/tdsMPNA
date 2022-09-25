#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

unsigned long long rdtsc(void);
int randxy(int x, int y);
double randreal();
void randvec(double* input, int size);
double mean(double *a, int size);
double std(double *a, double meanval, int size);
void printvec(int size, double* input);
void writemeasures (int nfunc, int nomeasures,double **nflop, double **meanvals, double **stdvals, char **fnames, char *file);

void saxpy(int size, double a, double* x, double* y)// 2
{
  #pragma omp for
	for(int i = 0; i<size; i++)
	{
		y[i] = a * y[i] + x[i];
	}
}

double dotprod (int size, double beta, double* x, double* y)
{
	double output = beta;
  #pragma omp reduction(+:output)
	for(int i = 0; i<size; i++)
	{
		output += y[i] * x[i];
	}
	return output;
}

double red (int size, double* x)
{
  double output = x[0];
  #pragma omp reduction(+:output)
  for(int i = 1; i<size; i++)
  {
    output +=x[i];
  }
	return	output;
}

int main(int argc, char** argv)
{
  if (argc != 4)
  {
    printf("invalid number of arguments ./prog maxsize numofmeasures filename\n");
    exit(0);
  }
  int maxsize = atoi(argv[1]), nomeasures = atoi(argv[2]), step =maxsize/nomeasures, sizemeasures = maxsize/step+1;
  int nrep = 10, i = 0, j = 0, k = 0, nfunc = 3;
  char * filename = argv[3];
  unsigned long long start = 0, end = 0, nocycles = 0;
  double meanvals[nfunc], stdvals[nfunc], nflop[nfunc], memory[nfunc];
  double ticks[nfunc][nrep];
  double *x, *y;
  double beta = randreal(), alpha = randreal(), rate = 0.0, res = 0.0;

  //file opening
  char fnames[3][10] = {"saxpy", "dotprod", "red(+)"} ;
  FILE *fp;
  fp = fopen(filename, "w+");
  fprintf(fp,"\n");
  for(int i = 0; i<nfunc; i++)
  {
      fprintf(fp,"%s; ; ; ; ;", fnames[i]);
  }
  fprintf(fp,"\n");
  for(int i = 0; i<nfunc; i++)
  {
      fprintf(fp,"no.elts;M.N.O.T;T.std; nflop;size(kB);");
  }
  fprintf(fp,"\n");

  for(i = step; i<=maxsize; i+=step)
  {
    x = (double*)calloc(i, sizeof(double));
    y = (double*)calloc(i, sizeof(double));
    nflop[0] = (double)(3*i);
    nflop[1] = (double)(3*(i-1));
    nflop[2] = (double)(i-1);
    memory[0] = ((double)(sizeof(double)*(i+1)))/1024.0;
    memory[1] = ((double)(sizeof(double)*(i+1)))/1024.0;
    memory[2] = ((double)(sizeof(double)*i))/1024.0;
    for(j = 0; j<nrep; j++)
    {
      randvec(x,i);
      randvec(y,i);
      #pragma omp parallel
      {
        #pragma omp single
        {
          start = rdtsc();
        }
        saxpy(i,alpha, x, y);
        #pragma omp single
        {
          end = rdtsc();
          printf("%d\n",k);
          printf("mean(y=ax+beta*y) = %lf\n", mean(y,i));
          ticks[0][j] = (double)(end -start);
          start = rdtsc();
        }
        res = dotprod (i, beta, x, y);
        #pragma omp single
        {
          end = rdtsc();
          ticks[1][j] = (double)(end -start);
          printf("dotprod = %lf\n",res);
          start = rdtsc();
        }
        res = red(i, y);
        #pragma omp single
        {
          end = rdtsc();
          ticks[2][j] = (double)(end -start);
          printf("red = %lf\n",res);

        }
      }
    }
    // saxpy
    meanvals[0] = mean(*(ticks), nrep);
    stdvals[0] = std(*(ticks), meanvals[0], nrep);
    //dotprod
    meanvals[1] = mean(*(ticks+1), nrep);
    stdvals[1] = std(*(ticks+1), meanvals[1], nrep);
    //reduction (/)
    meanvals[2] = mean(*(ticks+2), nrep);
    stdvals[2] = std(*(ticks+2), meanvals[2], nrep);
    for(int j = 0; j<nfunc; j++)
    {
        fprintf(fp,"%d;%f;%f;%f;%f;", i, meanvals[j], stdvals[j], nflop[j], memory[j]);
    }
    fprintf(fp,"\n");
    k++;
  }
  fclose(fp);
  free(x);
  free(y);


	return 0;
}

unsigned long long rdtsc(void)
{
  unsigned long long a, d;

  __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));

  return (d << 32) | a;
}

int randxy(int x, int y)
{
  return (rand() % (y - x + 1)) + x;
}

//
double randreal()
{
  int s = (randxy(0, 1)) ? 1 : -1;
  int a = randxy(1, RAND_MAX), b = randxy(1, RAND_MAX);

  return s * ((double)a / (double)b);
}

void randvec(double* input, int size)
{
  for(int i = 0; i<size; i++)
  {
    *(input+i) = randreal();
  }
}

double mean(double *a, int size)
{
	double r = 0.0;
	for(int i = 0; i < size; i++)
  {
    r += a[i];
  }
	return r / (double) size;
}

double std(double *a, double meanval, int size)
{
	double r = 0.0;
	for(unsigned int i = 0; i < size; i++)
  {
    r += a[i] * a[i];
  }
	r /= (double) size;
	r -= meanval*meanval;
	return sqrt(r);
}

void printvec(int size, double* input)
{
	printf("vector display (size = %d)\n", size);

	for(int i = 0; i<size; i++)
	{
		printf("vector[%d] = %lf\n",i,*(input+i));
	}
}

void printmatrix(int nolines, int nocols, double **input)
{
	printf("matrix display (%d \u2A09 %d)\n", size);

	for(int i = 0; i<size; i++)
	{
		printf("vector[%d] = %lf\n",i,*(input+i));
	}
}
