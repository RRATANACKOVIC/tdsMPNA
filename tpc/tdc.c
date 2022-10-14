#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

struct gso{//gso = gram schmidt output
  int n, m;
  double **Vm, **Hm;
};

unsigned long long rdtsc(void);
unsigned long long mean(unsigned long long *a, int size);
unsigned long long std(unsigned long long *a, unsigned long long meanval, int size);
double **initmat(int nolines, int nocols);
double *initvec(int length);
void freemat(double **mat, int nolines, int nocols);
void printmat(double **mat, int nolines, int nocols);
void printvec(double *vec, int length);
double dotprod(double *x, double *y, int length);
double norme_euclide(double *x, int length);
struct gso mgs(double **A, double *v, int n, int m);

int main (int argc, char **argv)
{
  if (argc != 6)
  {
    printf("wrong number of arguments: ./prog n m nomes funcname filename\n");
    exit(1);
  }

  int nrep = 10;
  int n = atoi(argv[1]), m = atoi(argv[2]), nomes = atoi(argv[3]), step = n/nomes;
  char *funcname = argv[4], *filename = argv[5];
  unsigned long long start = 0, end = 0, nocycles = 0;
  struct timespec tstart={0,0}, tend={0,0};

  unsigned long long meanval, stdval, nflop, memory;
  unsigned long long ticks[nrep];

  double **input, **outputgs;
  double *inputvec;
  struct gso outputcgs;

  FILE *fp;
  fp = fopen(filename, "w+");
  fprintf(fp,"\n");
  fprintf(fp,"%s; ; ; ; ;\n", funcname);
  fprintf(fp,"no.elts;M.N.O.T;T.std ; nflop;size(kB);\n");
  for(int i = step; i<=n; i+=step)
  {
    for(int j = 0; j<nrep; j++)
    {
      input = initmat(i,i);
      inputvec = initvec(i);
      start = rdtsc();
      outputcgs = mgs(input, inputvec, i, i);
      end = rdtsc();
    }

      //printmat(outputcgs.Hm, i, i+1);
      nflop = 2*(i*(i*i+1)+i*i*(i+1)+i*i);
      memory = 64*3*i*i+2*i;
    }
    fprintf(fp,"%d;%lld;%lld;%lld;%lld;\n", i, meanval, stdval, nflop, memory);
  }
  fclose(fp);

  return 0;
}

unsigned long long rdtsc(void)
{
  unsigned long long a, d;

  __asm__ volatile ("rdtsc" : "=a" (a), "=d" (d));

  return (d << 32) | a;
}

unsigned long long mean(unsigned long long *a, int size)
{
	unsigned long long r = 0;
	for(int i = 0; i < size; i++)
  {
    r += a[i];
  }
  r/=size;
	return r;
}

unsigned long long std(unsigned long long *a, unsigned long long meanval, int size)
{
	unsigned long long r = 0;
	for(int i = 0; i < size; i++)
  {
    r += a[i] * a[i];
  }
	r /= (unsigned long long) size;
	r -= meanval*meanval;
	return sqrt(r);
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


double **initmat(int nolines, int nocols)
{
  double **mat = (double**)malloc(nolines*sizeof(double*));
  for(int i = 0; i<nolines; i++)
  {
    mat[i] = (double*)calloc(nocols, sizeof(double));
  }
  for(int l = 0; l<nolines; l++)
  {
    for(int c = 0; c<nolines; c++)
    {
      mat[l][c] = randreal();
    }
  }
  return mat;
}

double *initvec(int length)
{
  double *output = calloc(length, sizeof(double));
  for(int i = 0; i<length; i++)
  {
    output[i] =randreal();
  }
  return output;
}

void freemat(double **mat, int nolines, int nocols)
{
  for(int i = 0; i<nolines; i++)
  {
    free(mat[i]);
  }
  free(mat);
}

void printmat(double **mat, int nolines, int nocols)
{
  for(int l = 0; l<nolines; l++)
  {
    for(int c = 0; c<nocols; c++)
    {
      printf("%lf ", mat[l][c]);
    }
    printf("\n");
  }
}

void printvec(double *vec, int length)
{
  for(int i = 0; i<length; i++)
  {
    printf("%lf ", vec[i]);
  }
  printf("\n");
}

double dotprod(double *x, double *y, int length)
{
  double output = 0.0;
  for (int i = 0; i<length; i++)
  {
    output+= x[i]*y[i];
  }
  return output;
}

double norme_euclide(double *x, int length)
{
  return sqrt(dotprod(x,x,length));
}


double *dpmv (double **A, double *x, int nolines, int nocols)
{
  double *output = (double*)calloc(nocols, sizeof(double));
  for(int i = 0; i<nolines; i++)
  {
    for(int j = 0; j<nocols; j++)
    {
      output[i] += x[j]*A[i][j];
    }
  }
  return output;
}


struct gso mgs(double **A, double *v, int n, int m)
{
  struct gso output;
  output.m = m;
  output.n = n;
  output.Hm = initmat(m,m+1);
  output.Vm = initmat(m+1,n);
  //double **Q = initmat(m+1,n);
  double * w;
  double b = 0.0, eps = 0.0000001, normpmo = 1/norme_euclide(v, n);//n
  for(int x = 0; x<n; x++)
  {
    output.Vm[0][x] = v[x]*normpmo;//n ->2n
  }
  for(int k = 1; k<m+1; k++)
  {
    w = dpmv(A,output.Vm[k-1],n,n);//2mn² -> 2n(mn+1)
    //printf("1\n");
    for(int j = 0; j<k; j++)
    {
      output.Hm[k-1][j] = dotprod(output.Vm[j], w, n);// 2n*m(m+1)/2 -> 2n(mn+1) +nm(m+1)
      for(int y = 0; y<n; y++)
      {
        w[y] -= output.Hm[k-1][j] *output.Vm[j][y];// 2n*m(m+1)/2 -> 2(n(mn+1)+nm(m+1))
      }
    }
    output.Hm[k-1][k] = norme_euclide(w,n);//2nm -> 2(n(mn+1)+nm(m+1)+mn)
    //printf("2\n");
    if(output.Hm[k-1][k]>eps)
    {
      for(int z = 0; z<n; z++)
      {
        output.Vm[k][z] = w[z]/output.Hm[k-1][k];
      }
    }
    else
    {
      return output;
    }
  }
  return output;
}
