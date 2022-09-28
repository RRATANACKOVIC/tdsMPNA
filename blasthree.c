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


void dgemm (struct mat *A, struct mat *B, struct mat *C, double alpha, double beta)
{
  // mem 64(i*j+3i+2+1)
  #pragma omp single
  {
    if( (A->nocols != B->nolines) || (A->nolines != C->nolines) || (B->nocols != C->nocols) )
    {
      printf("Dimension mismatch \n");
      exit(1);
    }
  }
  double temp = 0.0;
  struct mat TB;
  transmat(B,&TB);
  #pragma omp for
  for(int i = 0; i<A->nolines; i++)
  {
    for(int j = 0; j<B->nolines; j++)
    {
      temp = 0;
      for(int k = 0; k< A->nocols; k++)
      {
        temp += A->data[i][k]*TB->data[j][k];//2i
      }
      C->data[i][j] = alpha*temp + beta* C->data[i][j];// 3 i**2
    }
  }
}



void dsyrk (struct mat *A, struct mat *B, struct mat *C, double alpha, double beta)
{
  // mem 64(i*j+3i+2+1)
  #pragma omp single
  {
    if( (A->nocols != B->nolines) || (A->nolines != C->nolines) || (B->nocols != C->nocols) )
    {
      printf("Dimension mismatch \n");
      exit(1);
    }
  }
  struct mat TA, TB;
  //transmat(A, &TA);
  //transmat(B, &TB);
  #pragma omp for
  for(int i = 0; i<A->nolines; i++)
  {
    for(int j = 0; j<B->nolines; j++)
    {
      for(int k = 0; k< A->nocols; k++)
      {
        C->data[i][j] += alpha*A->data[i][k]*B.data[j][k] + beta*B->data[i][k]*A.data[j][k];//6 i**3
      }
    }
  }
}

double trace (struct mat *input)
{
  double output = 0.0;
  for (int i = 0; i<input->nolines; i++)
  {
    output += input->data[i][i];
  }
  return output;
}

int main(int argc, char** argv)
{
  if (argc != 4)
  {
    printf("invalid number of arguments ./prog maxsize numofmeasures filename\n");
    exit(0);
  }

  int maxsize = atoi(argv[1]), nomeasures = atoi(argv[2]), step =maxsize/nomeasures, sizemeasures = maxsize/step+1;
  int nrep = 10, i = 0, j = 0, k = 0, nfunc = 2;
  char * filename = argv[3];

  unsigned long long start = 0, end = 0, nocycles = 0;
  struct timespec tstart={0,0}, tend={0,0};

  double meanvals[nfunc], stdvals[nfunc], nflop[nfunc], memory[nfunc];
  double ticks[nfunc][nrep];
  double *x, *y;
  double beta = randreal(), alpha = randreal(), rate = 0.0, res = 0.0;
  struct mat A, TA, B, TB, C;

  //file opening
  char fnames[2][10] = {"dgemm", "dsyrk"};
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
      //fprintf(fp,"no.elts;M.N.O.T;T.std; nflop;size(kB);");
      fprintf(fp,"no.elts;M.RT(ns);RT.std (ns); nflop;size(kB);");
  }
  fprintf(fp,"\n");

  for(i = step; i<=maxsize; i+=step)
  {
    initmat(&A, i, i);
    initmat(&B, i, i);
    initmat(&C, i, i);
    nflop[0] = (double)(6*i*i+2);
    nflop[1] = nflop[0];
    memory[0] = ((double)(sizeof(double)*(2*i*i+2)))/1024.0;
    memory[1] = ((double)(sizeof(double)*(4*i*i+2)))/1024.0;
    for(j = 0; j<nrep; j++)
    {
      randmat(&A);
      randmat(&B);
      randmat(&C);
      #pragma omp parallel
      {
        #pragma omp single
        {
          //start = rdtsc();
          clock_gettime(CLOCK_MONOTONIC_RAW, &tstart);
        }
        dgemm (&A, &B, &C, alpha, beta);
        #pragma omp single
        {
          //end = rdtsc();
          clock_gettime(CLOCK_MONOTONIC_RAW, &tend);
          printf("%d\n",k);
          printf("dgemv = %lf\n", trace(&C));
          //ticks[0][j] = (double)(end -start);
          ticks[0][j] = (double)(tend.tv_nsec -tstart.tv_nsec);
          //start = rdtsc();
          clock_gettime(CLOCK_MONOTONIC_RAW, &tstart);
        }
        dgemm (&A, &B, &C, alpha, beta);
        #pragma omp single
        {
          //end = rdtsc();
          clock_gettime(CLOCK_MONOTONIC_RAW, &tend);
          //ticks[1][j] = (double)(end -start);
          ticks[1][j] = (double)(tend.tv_nsec -tstart.tv_nsec);
          printf("dsyrk = %lf\n",trace(&C));
        }
      }
    }
    // saxpy
    meanvals[0] = mean(*(ticks), nrep);
    stdvals[0] = std(*(ticks), meanvals[0], nrep);
    //dotprod
    meanvals[1] = mean(*(ticks+1), nrep);
    stdvals[1] = std(*(ticks+1), meanvals[1], nrep);
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
  freemat(&A);


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
  initmat(output, input->nolines, input->nocols);
  #pragma omp for
  for (int lin = 0; lin<input->nolines; lin++)
  {
    for (int cin = 0; cin<input->nocols; cin++)
    {
      output->data[cin][lin] = input->data[lin][cin];
    }
  }
}
