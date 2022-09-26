#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <time.h>

struct mat
{
  int nolines, nocols;
  double ** data;
};

unsigned long long rdtsc(void);
int randxy(int x, int y);
double randreal();
void randvec(double* input, int size);
double mean(double *a, int size);
double std(double *a, double meanval, int size);
void printvec(int size, double* input);
void writemeasures (int nfunc, int nomeasures,double **nflop, double **meanvals, double **stdvals, char **fnames, char *file);
void initmat (struct mat *input, int cols, int lines);
void printmat (struct mat *input);
void randmat (struct mat *input);
void freemat (struct mat *input);
void transmat (struct mat *input, struct mat * output);


void dgemv (struct mat *A, double *x, double *y, int sx, int sy, double alpha, double beta)
{
  // mem 64(i*j+3i+2+1)
  if(A->nocols != sx)
  {
    printf("Dimension mismatch (nocols of A and size of x)\n");
    exit(1);
  }
  if(A->nolines != sy)
  {
    printf("Dimension mismatch (nolines of A and size of y)\n");
    exit(1);
  }
  double z[sy];
  for(int i = 0; i<A->nolines; i++)
  {
    z[i] = 0.0;
    for(int j = 0; j<A->nocols; j++)
    {
      z[i] += x[j]*A->data[i][j];// 2*i*j
    }
    y[i] = alpha*z[i]+beta*y[i];// 3*j
  }
}

void dger (struct mat *A, double *x, double *y, int sx, int sy, double alpha)
{
  // mem 64(i*j+2i+1 +2)
  if(A->nolines != sx)
  {
    printf("Dimension mismatch (nocols of A and size of x)\n");
    exit(1);
  }
  if(A->nocols != sy)
  {
    printf("Dimension mismatch (nolines of A and size of y)\n");
    exit(1);
  }
  double z[sx];
  for(int k = 0; k<sx; k++)
  {
    z[k] = alpha *x[k];//i
    double z[sy];
  }
  for(int i = 0; i<A->nolines; i++)
  {
    for(int j = 0; j<A->nocols; j++)
    {
      A->data[i][j] = z[i]*y[j]+A->data[i][j];//2ij
    }
  }
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
  char fnames[3][10] = {"dgemv", "tdgemv", "dger"} ;
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
    nflop[0] = (double)(2*i*i+3*i);
    nflop[1] = nflop[0];
    nflop[2] = (double)(2*i*i+i);
    memory[0] = ((double)(sizeof(double)*(i*i+3*i+3)))/1024.0;
    memory[1] = ((double)(sizeof(double)*(2*i*i+3*i+3)))/1024.0;
    memory[2] = ((double)(sizeof(double)*(i*i+2*i+3)))/1024.0;
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
*/

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

void initmat (struct mat *input, int cols, int lines)
{
  //data_soa struct
  input->nolines = lines;
  input->nocols = cols;
  input->data = calloc(input->nolines, sizeof(double*));
  for (int l = 0; l<input->nolines; l++)
  {
    input->data[l] = calloc(input->nocols,sizeof(double));
  }
}

void freemat (struct mat *input)
{
  input->data = calloc(input->nolines, sizeof(double*));
  for (int l = 0; l<input->nolines; l++)
  {
    free(input->data[l]) ;
  }
  free(input->data);

}

void printmat (struct mat *input)
{
  for (int l = 0; l<input->nolines; l++)
  {
    for (int c = 0; c<input->nocols; c++)
    {
      printf("%lf ", input->data[l][c]);
    }
    printf("\n");
  }
}

void randmat (struct mat *input)
{
  for (int l = 0; l<input->nolines; l++)
  {
    for (int c = 0; c<input->nocols; c++)
    {
      input->data[l][c] = randreal();
    }
  }
}

void transmat (struct mat *input, struct mat * output)
{
  #pragma omp single
  {
      initmat(output, input->nolines, input->nocols);
  }
  #pragma omp for
  for (int lin = 0; lin<input->nolines; lin++)
  {
    for (int cin = 0; cin<input->nocols; cin++)
    {
      output->data[cin][lin] = input->data[lin][cin];
    }
  }
}
