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
double norme_frobenius(double **A, int nolines, int nocols);
double **transmat (double **input, int nolines, int nocols);
double *dpmv (double **A, double *x, int nolines, int nocols);
double **arnoldi(double **A, int n, int m);
struct gso mgs(double **A, double *v, int n, int m);
struct gso cgs(double **A, double *v, int n, int m);
double test_gs(double **Q, int nolines, int nocols);
double **LUPDecompose(double **A, int N, double Tol, int *P);
double LUPDeterminant(double **A, int *P, int N);
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
      if(strcmp(funcname,"gs") == 0)
      {
        input = initmat(i,i);
        start = rdtsc();
        outputgs = arnoldi(input, i, i);
        end = rdtsc();
        printf("sum(qi.qj) = %lf\n", test_gs(outputgs,i,i));
      }

      else if(strcmp(funcname,"mgs") == 0)
      {
        input = initmat(i,i);
        inputvec = initvec(i);
        start = rdtsc();
        outputcgs = mgs(input, inputvec, i, i);
        end = rdtsc();
        printf("sum(qi.qj) = %lf\n", test_gs(outputcgs.Vm,i+1,i));
      }
      else if(strcmp(funcname,"cgs") == 0)
      {
        input = initmat(i,i);
        inputvec = initvec(i);
        start = rdtsc();
        outputcgs = cgs(input, inputvec, i, i);
        end = rdtsc();
        printf("size, %d, sum(qi.qj) = %lf\n",i, test_gs(outputcgs.Vm,i+1,i));
      }
      ticks[j] = end - start;
    }
    meanval = mean(ticks, nrep);
    stdval = std(ticks, meanval, nrep);
    if(strcmp(funcname,"gs") == 0)
    {
      memory = 64*2*i*(i+1) + 2;
      nflop = i*(i+1)*(4*i+1)/2 + i*i + 2*i + 1;
    }

    else if(strcmp(funcname,"mgs") == 0)
    {
      //printmat(outputcgs.Hm, i, i+1);
      nflop = 2*(i*(i*i+1)+i*i*(i+1)+i*i);
      memory = 64*3*i*i+2*i;
    }

    else if(strcmp(funcname,"cgs") == 0)
    {
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


double **arnoldi(double **A, int n, int m)
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


struct gso cgs(double **A, double *v, int n, int m)
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

    }
    for(int h = 0; h<k; h++)
    {
      for(int y = 0; y<n; y++)
      {
        w[y] -= output.Hm[k-1][h] *output.Vm[h][y];// 2n*m(m+1)/2 -> 2(n(mn+1)+nm(m+1))
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

double test_gs(double **Q, int nolines, int nocols)
{
  double output = 0.0;
  for(int i = 0; i<nolines; i++)
  {
    for(int j = i+1; j<nolines; j++)
    {
      output += dotprod(Q[i], Q[j], nocols);
    }
  }
  return output;
}

double **LUPDecompose(double **A, int N, double Tol, int *P)
{
    double **output;
    int i, j, k, imax;
    double maxA, *ptr, absA;

    for (i = 0; i <= N; i++)
        P[i] = i; //Unit permutation matrix, P[N] initialized with N

    for (i = 0; i < N; i++) {
        maxA = 0.0;
        imax = i;

        for (k = i; k < N; k++)
            if ((absA = fabs(A[k][i])) > maxA) {
                maxA = absA;
                imax = k;
            }

        if (maxA < Tol) return 0; //failure, matrix is degenerate

        if (imax != i)
        {
            //pivoting P
            j = P[i];
            P[i] = P[imax];
            P[imax] = j;

            //pivoting rows of A
            ptr = A[i];
            A[i] = A[imax];
            A[imax] = ptr;

            //counting pivots starting from N (for determinant)
            P[N]++;
        }

        for (j = i + 1; j < N; j++) {
            A[j][i] /= A[i][i];

            for (k = i + 1; k < N; k++)
                A[j][k] -= A[j][i] * A[i][k];
        }
    }

    return output;  //decomposition done
}

double LUPDeterminant(double **A, int *P, int N)
{

    double det = A[0][0];

    for (int i = 1; i < N; i++)
        det *= A[i][i];

    return (P[N] - N) % 2 == 0 ? det : -det;
}
