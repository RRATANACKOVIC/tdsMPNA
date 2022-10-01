#include <stdio.h>
#include <stdlib.h>

void initmat(double **mat, int nolines, int nocols);
void freemat(double **mat, int nolines, int nocols);
void printmat(double **mat, int nolines, int nocols);
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
  double **summat;
  initmat(summat, nolines, nocols);
  printf(" 1\n");
  printmat(summat, nolines, nocols);
  printf(" 2\n");
  freemat(summat, nolines, nocols);
  printf(" 3\n");
  return 0;
}


void initmat(double **mat, int nolines, int nocols)
{
  mat = (double**)malloc(nolines*sizeof(double*));
  for(int i = 0; i<nolines; i++)
  {
    mat[i] = (double*)calloc(nocols, sizeof(double));
  }
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
    for(int c = 0; c<nocols; c ++)
    {
      printf("%lf ", mat[l][c]);
      printf("%lf ", *(*(mat+l)+c));
    }
    printf("\n");
  }
}
