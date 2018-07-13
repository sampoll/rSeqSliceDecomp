#include <R.h>
#include <Rdefines.h>
#include <Rinternals.h>

#include "slicefinder.h"
#include <stdlib.h>
#include <string.h>

int check(unsigned *Y, slice **slices, unsigned ns, unsigned *rpoints, unsigned nr, unsigned np);

SEXP slicify(SEXP rx, SEXP rnp)  {

  int *np_p = INTEGER(rnp);
  int *xp_p = INTEGER(rx);

  unsigned *X = (unsigned *)malloc((*np_p)*sizeof(unsigned));
  unsigned np = (unsigned) (*np_p);
  for(int i=0;i<np;i++)  {
    X[i] = (unsigned) xp_p[i]; 
  }

  // make a copy of X to check afterwards
  unsigned *Y = (unsigned *)malloc(np*sizeof(unsigned));
  for(int i=0;i<np;i++)  {
    Y[i] = X[i];
  }

  // slice decomposition
  int nrpoints = -1, nslices = -1;
  unsigned *rpoints = NULL;
  slice **slices = NULL;

  algorithm1(X, np, &rpoints, &nrpoints, &slices, &nslices);

  if (!check(Y, slices, nslices, rpoints, nrpoints, np))  {
    printf("error: slice decomposition\n");
  }

  int ntotal = nslices + nrpoints / 2;
  if (nrpoints % 2)
    ntotal++;

  char **sres = (char **)malloc(ntotal*sizeof(char *));
  int j = 0;
  for(int i=0;i<nslices;i++)  {
    size_t sz = snprintf(NULL, 0, "%u:%u:%u", slices[i]->start, slices[i]->stop, slices[i]->step);
    sres[j] = (char *)malloc((sz+1)*sizeof(char));
    if (sres[j] == NULL)  {
      printf("error: out of memory!\n");
    }
    sprintf(sres[j], "%u:%u:%u", slices[i]->start, slices[i]->stop, slices[i]->step);
    j++;
  }
  for(int i=0;i<nrpoints/2;i++)  {
    unsigned start = rpoints[2*i];
    unsigned stop = rpoints[2*i+ 1];
    unsigned step = stop-start;

    size_t sz = snprintf(NULL, 0, "%u:%u:%u", start, stop, step);
    sres[j] = (char *)malloc((sz+1)*sizeof(char));
    if (sres[j] == NULL)  {
      printf("error: out of memory!\n");
    }
    sprintf(sres[j], "%u:%u:%u", start, stop, step);
    j++;
  }
  if (nrpoints % 2)  {
    unsigned start = rpoints[nrpoints-1];
    unsigned stop = start;
    unsigned step = 1;
    
    size_t sz = snprintf(NULL, 0, "%u:%u:%u", start, stop, step);
    sres[j] = (char *)malloc((sz+1)*sizeof(char));
    if (sres[j] == NULL)  {
      printf("error: out of memory!\n");
    }
    sprintf(sres[j], "%u:%u:%u", start, stop, step);
    j++;
  }

  for(int i=0;i<nslices;i++)  {
    free(slices[i]->S);
    free(slices[i]);
  }
  free(slices);
  free(rpoints);
  free(Y);

  SEXP rsres = Rf_allocVector(STRSXP, ntotal);
  PROTECT(rsres);
  for(int i=0;i<ntotal;i++)  {
    SET_STRING_ELT(rsres, i, mkChar(sres[i]));
  }
  UNPROTECT(1);
  return(rsres);

}




// Temporary: check that no points are lost and no spurious points gained

int compare(const void *a, const void *b)  {  return ( *(unsigned *)a - *(unsigned *)b );  }
int check(unsigned *Y, slice **slices, unsigned ns, unsigned *rpoints, unsigned nr, unsigned np)  {
  unsigned *Z = (unsigned *)malloc(np*sizeof(unsigned));
  int nz = 0;
  unsigned toomany = 0;

  for(int i=0;i<ns;i++)  {
    for(int s=slices[i]->start;s<=slices[i]->stop;s += slices[i]->step)  {
      if (nz < np)
        Z[nz++] = s;
      else
        toomany = 1;
    }
  }
  for(int i=0;i<nr;i++)  {
      if (nz < np)
        Z[nz++] = rpoints[i];
      else
        toomany = 1;
  }

  qsort(Z, nz, sizeof(unsigned), compare);

  int ok = 1;

  if (toomany)  {
    printf("error: too many points\n");
    ok = 0;
  }
  else if (nz != np) {
    printf("error: too few points\n");
    ok = 0;
  }
  else  {
    for(int i=0;i<np;i++)  {
      if (Y[i] != Z[i])
        ok = 0;
    }
  }

  free(Z);
  return(ok);

}


