#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct gso{//gso = gram schmidt output
  int n, m;
  double **Vm, **Hm;
};

double **initmat(int nolines, int nocols);
double *initvec(int length);
void freemat(double **mat, int nolines, int nocols);
void printmat(double **mat, int nolines, int nocols);
void printvec(double *vec, int length);
double dotprod(double *x, double *y, int length);
double norme_euclide(double *x, int length);
double norme_frobenius(double **A, int nolines, int nocols);
struct gso cgs(double **A, double *v, int n, int m);
int main (int argc, char **argv)
{
  if (argc != 3)
  {
    printf("wrong number of arguments: ./prog nolines nocols\n");
    exit(1);
  }

  int nolines = atoi(argv[1]), nocols = atoi(argv[2]);
  double **summat = initmat(nolines, nocols);
  double *sumvec = initvec(nocols);
  struct gso result = cgs(summat, sumvec, nolines, nocols);
  printmat(result.Hm, result.n, result.m);
  printmat(result.Vm, result.n, result.m);
  free(sumvec);
  freemat(summat, nolines, nocols);
  return 0;
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

double **getQ(double **A, double *v, int n, int m)
{
  double q = initmat(m,n);
  //for(int l = 0; )
}

double **gs(double **A, double *v, int n, int m)
{
  double **output = (double**)malloc(n*sizeof(double*));
  double dp = 0.0;
  for(int i = 0; i<n; i++)
  {
    output[i] = (double*)calloc(m, sizeof(double));
  }
  double normpmo = 1.0/norme_euclide(A[0], n);
  for(int i = 0; i<n; i++)
  {
    output[0][i] = normpmo * A[0][i];
  }
  for(int i = 1; i<m; i++)
  {
    for(int p = 0; p<n; p++)
    {
      output[i][p] = A[i][p];
    }
    for(int j = 0; j<j-1; j++)
    {
      dp = dotprod(output[i], output[i],  n);
      for(int q = 0; q<n; q++)
      {
        output[i][q] -= dp*output[j][q];
      }
    }
    normpmo = 1.0/norme_euclide(A[i], n);
    for(int t = 0; t<n; t++)
    {
      output[i][t] *= normpmo;
    }
  }

}

struct gso cgs(double **A, double *v, int n, int m)
{
  struct gso output;
  output.m = m;
  output.n = n;
  output.Hm = initmat(n,m);
  output.Vm = initmat(n,m);
  double *qj = initvec(n);
  double * w = initvec(n);
  for(int k = 0; k<n; k++)
  {
    for(int i = 0; i<n; i++)
    {
      w[i] = A[i][k];
    }
    for(int j = 0; j<k; k++)
    {
      output.Hm[j][k] = dotprod (w,qj,n);
    }
    for(int i = 0; i<n; i++)
    {
      for(int j = 0; j<k; k++)
      {
        w[i]-= output.Hm[j][k]*qj[i];
      }
    }
    output.Hm[k][k] = norme_euclide (w, n);
    for(int i = 0; i<n; i++)
    {
      qj[i] = w[i]/output.hm[k][k];
    }
  }

  return output;
}
