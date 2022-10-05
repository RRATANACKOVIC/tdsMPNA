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
double std(unsigned long long *a, unsigned long long meanval, int size);
double **initmat(int nolines, int nocols);
double *initvec(int length);
void freemat(double **mat, int nolines, int nocols);
void printmat(double **mat, int nolines, int nocols);
void printvec(double *vec, int length);
double dotprod(double *x, double *y, int length);
double norme_euclide(double *x, int length);
double norme_frobenius(double **A, int nolines, int nocols);
double **transmat (double **input, int nolines, int nocols);
double *dpmv (double **A, double *x, int nolines, int nocols);
double **gs(double **A, int n, int m);
struct gso cgs(double **A, double *v, int n, int m);
int main (int argc, char **argv)
{
  if (argc != 6)
  {
    printf("wrong number of arguments: ./prog n m nomes funcname filename\n");
    exit(1);
  }

  int nrep = 10;
  int n = atoi(argv[1]), m = atoi(argv[2]), nomes = atoi(argv[3]), step = n/nomes;
  char *funcname = argv[4];
  unsigned long long start = 0, end = 0, nocycles = 0;
  struct timespec tstart={0,0}, tend={0,0};

  unsigned long long meanval, stdval, nflop, memory;
  unsigned long long ticks[nrep];

  FILE *fp;
  fp = fopen(filename, "w+");
  fprintf(fp,"\n");
  fprintf(fp,"%s; ; ; ; ;\n", funcname);
  fprintf(fp,"no.elts;M.N.O.T;T.std ; nflop;size(kB);\n");
  for(int i = 0; i<n; i+=step)
  {
    for(int j = 0; j<nrep; j++)
    {
      start = rdtsc();
      //runfunc
      end = rdtsc();
      ticks[i] = end - start;
    }
    meanval = mean(ticks, nrep);
    stdval = std(ticks, meanval, nrep);
    memory = 2*i*(i+1) + 2;
    nflop = i*(i+1)*(4*i+1)/2 + i*i + 2*i + 1;
    fprintf(fp,"%d;%lld;%lf;%lld;%lld;\n", i, meanvals, stdvals, nflop, memory);
  }
  fclose(fp);
  /*
  double **summat = initmat(nolines, nocols);
  printmat(summat, nolines, nocols);
  printf("\n");
  double *sumvec = initvec(nocols);
  printvec(sumvec, nocols);
  printf("\n");
  double *davec = dpmv(summat, sumvec, nolines, nocols);
  printvec(davec, nolines);
  printf("\n");
  */
  /*
  double *sumvec = initvec(nocols);
  struct gso result = cgs(summat, sumvec, nolines, nocols);
  printmat(result.Hm, result.n, result.m);
  printmat(result.Vm, result.n, result.m);
  free(sumvec);
  */
  /*
  printmat(summat, nolines, nocols);
  double ** omat = gs(summat, nolines, nocols);
  printmat(omat, nolines, nocols);
  freemat(omat, nolines, nocols);
  freemat(summat, nolines, nocols);
  */
  /*
  double **S = (double **)malloc(2*sizeof(double*));
  S[0] = (double *)calloc(2, sizeof(double));
  S[1] = (double *)calloc(2, sizeof(double));
  S[0][0] = 3.0;
  S[0][1] = 2.0;
  S[1][0] = 1.0;
  S[1][1] = 2.0;
  double **GS = gs(S,2,2);
  printmat(GS,2,2);
  */
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
	return r / (unsigned long long) size;
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

double norme_frobenius(double **A, int nolines, int nocols)
{
  double output = 0.0;
  for(int i = 0; i<nolines; i++)
  {
    for(int j = 0; j<nocols; j++)
    {
      output+= A[i][j]*A[i][j];
    }
  }
}

double **transmat (double **input, int nolines, int nocols)
{
  double **output = initmat(nocols, nolines);
  for (int lin = 0; lin<nolines; lin++)
  {
    for (int cin = 0; cin<nocols; cin++)
    {
      output[cin][lin] = input[lin][cin];
    }
  }
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


double **gs(double **A, int n, int m)
{
  // A in R(n*m) -> output ~ transpos(A) in R(m*n)
  double **output = (double**)malloc(m*sizeof(double*));
  double dp = 0.0, normpmo = 0.0;
  for(int x = 0; x<m; x++)
  {
    output[x] = (double*)calloc(n, sizeof(double));
  }
  for(int y = 0; y<n; y++)
  {
    output[0][y] = A[y][0];
  }
  normpmo = 1.0/norme_euclide(output[0], n);//1
  for(int y = 0; y<n; y++)
  {
    output[0][y] *= normpmo;//n ->n+1
  }
  for(int i = 1; i<m; i++)
  {
    for(int p = 0; p<n; p++)
    {
      output[i][p] = A[p][i];
    }
    for(int j = 0; j<i; j++)
    {
      dp = dotprod(output[i], output[j],  n); //m(m+1)/2 *(2n-1) -> m(m+1)(2n+1)/2 + n + 1
      for(int q = 0; q<n; q++)
      {
        output[i][q] -= dp*output[j][q]; // m(m+1)/2 *2n -> m(m+1)(4n+1)/2 + n + 1
      }
    }
    normpmo = 1.0/norme_euclide(output[i], n);//m -> m(m+1)(4n+1)/2 + m + n + 1
    for(int t = 0; t<n; t++)
    {
      output[i][t] *= normpmo;// mn -> m(m+1)(4n+1)/2 + mn + m + n + 1
    }
  }
  return output;
}

struct gso cgs(double **A, double *v, int n, int m)
{
  struct gso output;
  output.m = m;
  output.n = n;
  output.Hm = initmat(m,m+1);
  output.Vm = initmat(n,m);
  double **Q = initmat(m+1,n);
  double * q = initvec(n);
  double * w;
  double b = 0.0, eps = 0.0000001, normpmo = 1/norme_euclide(v, n);
  for(int x = 0; x<n; x++)
  {
    v[x] *=normpmo;
  }
  for(int k = 1; k<n+1; k++)
  {
    w = dpmv(A,Q[k-1],n,n);
    for(int j = 0; j<k; j++)
    {
      h[j][k-1] = dotprod(Q[j], w, n);
      for(int y = 0; y<n; y++)
      {
        v[y] -= output.Hm[j][k-1] *Q[j][y];
      }
    }
    output.Hm[k][k-1] = norme_euclide(v,n);
    if(output.Hm[k][k-1]>eps)
    {
      for(int z = 0; z<n; z++)
      {
        Q[k][z] = v[z]/output.Hm[k][k-1];
      }
    }
    else
    {
      return output;
    }
  }

}
