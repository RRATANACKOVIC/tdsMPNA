#include <stdio.h>
#include <stdlib.h>
#include <math.h>
void initmat(double *mat, int nolines, int nocols);
void freemat(double *mat, int nolines, int nocols);
void printmat(double *mat, int nolines, int nocols);
double dotprod(double *x, double *y, int length);
double norme_frobenius(double *A, int nolines, int nocols);
int main (int argc, char **argv)
{
  /*
  if (argc != 3)
  {
    printf("wrong number of arguments: ./prog nolines nocols\n");
    exit(1);
  }
  */
  int nolines = atoi(argv[1]), nocols = atoi(argv[2]);
  printf("%d %d\n", nolines, nocols);
  double *summat = calloc(nolines*nocols, sizeof(double));;
  initmat(summat, nolines, nocols);
  printmat(summat, nolines, nocols);
  printf("%lf\n", norme_frobenius(summat, nolines, nocols));
  free(summat);
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

void initmat(double *mat, int nolines, int nocols)
{
  for(int i = 0; i<nolines; i++)
  {
    for(int j = 0; j<nocols; j++)
    {
      mat[i*nocols + j] = randreal();
    }
  }
}



void printmat(double *mat, int nolines, int nocols)
{
  for(int l = 0; l<nolines; l++)
  {
    for(int c = 0; c<nocols; c ++)
    {
      //printf("%lf ", mat[l*nocols +c]);
      printf("%lf ", *(mat +l*nocols+ c));
    }
    printf("\n");
  }
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

double norme_frobenius(double *A, int nolines, int nocols)
{
  double output = 0.0;
  for(int i = 0; i<nolines; i++)
  {
    for(int j = 0; j<nocols; j++)
    {
      output+= A[i*nocols + j]*A[i*nocols + j];
    }
  }
}
