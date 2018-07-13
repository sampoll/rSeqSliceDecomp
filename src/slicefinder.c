#include "slicefinder.h"

// private
int pathlen(unsigned **D, unsigned nd, unsigned i, unsigned j);
int printpth(unsigned **D, unsigned *X, unsigned nd, unsigned ii, unsigned jj);

// private: list struct for building slice list
typedef struct _slicenode {
  slice *slice;
  struct _slicenode *nxt;
} slicenode;


// insert slice at top of list 
void insert_slice_list(slice *slp, slicenode **headp)  {
  slicenode *node = (slicenode *)malloc(sizeof(slicenode));
  node->slice = slp;
  node->nxt = *headp;
  *headp = node;
}

// delete slice list nodes but not the slices
void delete_slice_list(slicenode *head)  {
  slicenode *node = head;
  while (node != NULL)  {
    if (node->nxt == NULL)  {
      free(node);
      node = NULL;
    }
    else  {
      slicenode *nxtnode = node->nxt;
      free(node);
      node = nxtnode;
    }
  }
}

// Moderately Hard: find all slices of a given step size 
//   look through D for the value and find slices 
//   beginning with those (i,j) pairs

// Easy: find all slices of a given start index (or value)
// Must pass NULL, 0 to (slice, nr) or could crash
int allslices_start_given(sequence *Q, slice ***R, unsigned *nr, unsigned start)  {
  unsigned nd = Q->nd;
  unsigned istart, j, k;

  istart = -1;
  for(int j=0;j<Q->np;j++)  {
    if (Q->X[j] == start)  {
      istart = j;
      break;
    }
  }

  if (istart == -1)  {
    printf("error: start value %u is not in the sequence\n", start);
    return 1;
  }
 
  slicenode *head = NULL;
  for(j=0;j<nd-istart;j++)  {
    unsigned pl = 1 + pathlen(Q->D, nd, istart, j);   // number of intervals
    // every pair is a trivial slice, omit these
    if (pl >= 2)  {     
      slice *slp = initslice(Q, istart, j, pl);
      insert_slice_list(slp, &head);
    }
  }

  // count slices in list
  slicenode *node = head;
  for(j = 0;node != NULL;j++)  {
    node = node->nxt;
  }

  // no nontrivial slices begin at 
  if (j == 0)  {
    // printf("info: no nontrivial slices begin at %u\n", start);
    return 0;
  }

  // allocate array and put slices in it
  *nr = j;
  slice **RR = (slice **)malloc(j*sizeof(slice *));
  node = head;
  // put slices in backwards for increasing step size
  for(k=0;k<j;k++)  {
    RR[j-1-k] = node->slice;
    node = node->nxt;
  }
  *R = RR;
  
  // clean up list
  delete_slice_list(head);
  return 0;

}

int remove_subsets(slice ***L, unsigned *I, unsigned np)  {
    unsigned nsubsets = 0;
    // remove subsets
    for(unsigned istart=0;istart<np;istart++)  {
      for(unsigned j=0;j<I[istart];j++)  {
        slice *slp = L[istart][j];
        if (slp != NULL)  {         // NULL indicates an already-removed subset

          // slp->S has the indices of the slice
          // so subset slices are at S[1], S[2], ..., S[ns-3]
          // because trivial (i.e. two-point) slices are not
          // included in the lists

          for(unsigned k=1;k<slp->ns-2;k++)  {
            // find slice in list L[k]
            for(unsigned l=0;l<I[slp->S[k]];l++)  {
              slice *slq = L[slp->S[k]][l];
              if (slq != NULL && 
                  slq->start == slp->start + k*slp->step &&
                  slq->step == slp->step)  {
                // remove slice
                free(slq->S);
                free(slq);
                L[slp->S[k]][l] = NULL;
                nsubsets++;
                break;   // only one slice per (step, start) pair
              }
            }

          }
        }
      }
    }
    return(nsubsets);
}

int allslices(sequence *Q, slice ***C, unsigned *nc, unsigned include_subsets)  {

  // build one list for every start position, to facilitate removing 
  // subset slices. If slice S is start:stop:step, then its subsets are 
  // (start+step):stop:step, (start+2*step):stop:step, etc.
  unsigned *X = Q->X;
  int np = Q->np;
  int irtn = 0;

  slice ***L = (slice ***)malloc(np*sizeof(slice **));
  unsigned *I = (unsigned *)malloc(np*sizeof(unsigned));

  slice **R = NULL;
  unsigned nr;

  for(unsigned istart=0;istart<np;istart++)  {
    unsigned start = X[istart];
    
    R = NULL;
    nr = 0;
    irtn = allslices_start_given(Q, &R, &nr, start);

    L[istart] = R;
    I[istart] = nr;

  }

  unsigned nb = 0;
  if (!include_subsets)  {
    nb = remove_subsets(L, I, np);
  }

  // pack all slice arrays into one array
  unsigned sz = 0;
  for(unsigned istart=0;istart<np;istart++)  {
    sz += I[istart];
  }
  sz -= nb;    // subtract off space of removed subset slices
  *C = (slice **)malloc(sz*sizeof(slice *));

  unsigned i = 0;
  for(unsigned istart=0;istart<np;istart++)  {
    for(unsigned j=0;j<I[istart];j++)  {
      if (L[istart][j] != NULL)  {
        (*C)[i++] = L[istart][j];
      }
    }
  }
  *nc = sz;
  
  // delete lists in L but not the slices in them
  for(unsigned istart=0;istart<np;istart++)  {
    free(L[istart]);
  }
  free(L);
  free(I);

  return 0;
}



// Initial computation of all diffs from sequence points
sequence *initseq(unsigned *X, unsigned np)  {
  sequence *Q = (sequence *)malloc(sizeof(sequence));
  Q->X = X;
  Q->np = np;
  Q->nd = np-1;
  int i, j;

  int nd = Q->nd;    // convenience

  // D[i][j] is distance from X[i] to X[i+(j+1)]
  unsigned **D = (unsigned **)malloc(nd*sizeof(unsigned *));
  for(i=0;i<nd;i++)  {
    D[i] = (unsigned *)malloc(nd*sizeof(unsigned));     // really only need i
  }

  // compute diffs
  for(i=0;i<nd;i++)  
    D[i][0] = X[i+1] - X[i];

  for(i=1;i<nd;i++)  {
    for(j=0;j<nd-i;j++)  
      D[j][i] = D[j+1][i-1] + D[j][0];
    for(j=nd-i;j<nd;j++) 
      D[j][i] = 0;
  }

  Q->D = D;
  Q->S = NULL;
  return Q;
}

void printseq(sequence *Q)  {

  unsigned **D = Q->D;
  unsigned nd = Q->nd;
  unsigned np = Q->np;

  printf("Sequence of length %u\n\n", np);

  // view sequence
  for(int i=0;i<np;i++)  {
    printf("%u ", Q->X[i]);
  }
  printf("\n\n");

  // view diffs
  for(int i=0;i<nd;i++)  {
    for(int j=0;j<nd-i;j++)  {
      printf("%10u", D[i][j]);
    }
    printf("\n");
  }

}


// Read sequence from file. Eventually, initialize sequence from a SEXP.
int input(const char *fn, unsigned **xp, int *n)  {

  FILE *fp = fopen(fn, "r");
  if (fp == NULL)  {
    printf("error: can't open file %s\n", fn);
    exit(1);
  }

  unsigned nn = 0;
  int nnn;
  fscanf(fp, "%d", &nnn);

  nn = (unsigned) nnn;
  unsigned *X = (unsigned *)malloc(nn*sizeof(unsigned));

  for(int i=0;i<nn;i++)  {
    fscanf(fp, "%d", &nnn);
    X[i] = (unsigned) nnn;
  }
  fclose(fp);
  *n = nn;
  *xp = X;

  return 0;

}

// The sequence is easy to compute from the indices, but the 
// indices are not easy to compute from the points.
// I.e., given S = {i1, i2, i3, ... } it's easy to get the 
// slice points = {X[S[i1]], X[S[i2]], ... }
// But given points {X1, X2, ... } it's necessary
// to do a search to find i1, i2, etc. 

// Return length of slice from X[i] with step D[i][j] 
// D[i][j] is the difference between X[i] and X[i+j+1]

int pathlen(unsigned **D, unsigned nd, unsigned i, unsigned j)  {
  unsigned val = D[i][j];

  if ((i + j + 1) > nd) 
    return 0;

  // special case: last point of sequence is in slice
  if ((i + j + 1) == nd)  {   
    // return 1;
    return 0;
  }

  // if the step size (val) is in row (i+j+1)
  // then the length of the slice is one more
  // then the length of the rest of the slice.

  for(int k=0;k<(nd - (i+j+1));k++)  {
    if (D[i+j+1][k] == val)  {
      return (1 + pathlen(D, nd, i+j+1, k));
    }
  }
  return 0;

}

// private function for debugging only
int printpth(unsigned **D, unsigned *X, unsigned nd, unsigned ii, unsigned jj)  {

  unsigned i = ii, j = jj;
  unsigned val = D[i][j];
  unsigned found = 0;
  printf("%4u ", X[i]);
  i = (i + j + 1);

  while (1)  {
    if (i > nd)  {
      printf("warning: this should not happen\n");
      return 0;
    }

    if (i == nd)  { 
      printf("%4u ", X[nd]);
    }

    // find column number for diff val in row i
    found = 0;
    for(int k=0;found == 0 && k<nd-i;k++)  {
      if (D[i][k] == val)  {
        printf("%4u ", X[i]);
        i = i + k + 1;
        found = 1;
      }
    }

    if (!found)  {
      printf("%4u ", X[i]);
      break;
    }

  }
  printf("\n");
  return 0;

}

// Find the slice {start = X[i], step = (X[i+j]-X[i])}
// Need to use pathlen first to find the length to allocate
slice *initslice(sequence *Q, unsigned ii, unsigned jj, unsigned len)  {

  unsigned **D = Q->D;
  int nd = Q->nd;
  unsigned step = D[ii][jj];
  unsigned found = 0;
  unsigned i = ii, j = jj;

  // len is the number of steps; number of points is one more than that.
  slice *slc = (slice *)malloc(sizeof(slice));
  unsigned *S = (unsigned *)malloc((len+1)*sizeof(unsigned));

  S[0] = i;
  int s = 1;
  i = (i + j + 1);

  while (1)  {
    if (i > nd)  { 
      printf("warning: this should not happen\n");
      return 0;
    }

    // special case: last point in slice
    // if (i == nd)  {
    //  S[s++] = nd;
    // }

    found = 0;
    for(int k=0;found == 0 && k<nd-i;k++)  {
      if (D[i][k] == step)  {
        S[s++] = i;
        i = i + k + 1;
        found = 1;
      }
    }

    if (!found)  {
      S[s++] = i;
      break;
    }
  }

  slc->start = Q->X[S[0]];
  slc->step = step;
  slc->stop = Q->X[S[s-1]];
  slc->ns = s;
  slc->S = S;

  if (s != len+1)  {
    printf("error: found %u points but len is %u:  ", s, len+1);
    printf("start:stop:step = %u:%u:%u\n", slc->start, slc->stop, slc->step);
  }

  return slc;
}

void printslice(sequence *Q, slice *slc)  {
  printf("start:stop:step = %u:%u:%u\n", slc->start, slc->stop, slc->step);
  for(int i=0;i<slc->ns;i++)  
    printf("%u ", slc->S[i]);
  printf("  [");
  for(int i=0;i<slc->ns;i++)  
    printf("%u ", Q->X[slc->S[i]]);
  printf("]\n\n");
}


// Remove slice from Q 

// TODO: need to remove points from Q->X and update Q->np as well as Q->nd 
// NOTE: that when removing a slice, we remove all subset slices simultaneously,
//       Subsets need to be removed from array of slices before calling this.

// Assumptions: slice is at least 3 points (non-trivial slice)

int removeslice(sequence *Q, slice *slc)  {  
  unsigned **D = Q->D;
  unsigned *S = slc->S;
  unsigned nd = Q->nd;
  unsigned ns = slc->ns;

  // Remove relevant elements from difference array D
  // Removing point X[p] requires:
  // (2) Remove row p altogether
  // (1) Diagonally: remove all differences to X[p]
  //     I.e., element p in row 0, element p-1 in row 1, ...
  //     one element in each row up to p-1
  // Set elements to be removed to 0.

  if (Q->np == ns)  {  // done
    free(Q->X); 
    Q->X = NULL;
    Q->nd = 0;
    Q->np = 0;
    return 0;
  }

  // remove diagonally
  for(int i=0;i<ns;i++)  {
    int s = S[i];
    for(int l=0;l<s;l++)  { 
      D[l][s-l-1] = 0;
    }
  }

  // squeeze out obsolete rows
  // instead of copying elements just copy pointers. 
  int m = 0;    // pull up this many rows
  int j = 0;    // index into S
  for(int i=0;i<nd;i++)  {
    if (j < ns && S[j] == i)  {
      j++;
      m++;
      free(D[i]);
    }
    else  {
      D[i-m] = D[i];
    }
  }

  // squeeze out zeroes in remaining rows
  for(int i=0;i<nd-ns;i++)  {
    int m = 0;
    for(int k=0;k<nd-ns-i;k++)  {
      if (D[i][k] == 0)  {
        while(D[i][k+m] == 0)
          m++;
        D[i][k] = D[i][k+m];
        if (m > 0)
          D[i][k+m] = 0;
      }
    }
  }

  // TODO: remove this check eventually 
  // debug: check for non-zeroed out elements, because the 
  // zero is used as a signal of no-data

  for(int i=0;i<nd-ns;i++)  {
    for(int j=nd-ns;j<nd;j++)  {
      if (D[i][j] != 0)  
        printf("Non-zero (%u,%u)\n", i, j);
    }
  }

  // Remove slice points from Q->X
  unsigned *XX = (unsigned *)malloc((Q->np-ns)*sizeof(unsigned));
  unsigned ii = 0;
  unsigned jj = 0;
  for(int i=0;i<Q->np;i++)  { 
    if (ii < ns && i == S[ii])   // skip this point
      ii++; 
    else
      XX[jj++] = Q->X[i];
  }
  free(Q->X); 
  Q->X = XX;
  Q->nd = nd - ns;
  Q->np = Q->np - ns;
  return 0;

}

// private for algorithm1
int maxslice(slice **C, unsigned nc)  {
  unsigned maxns = 0;
  int maxidx = -1;

  if (nc == 0)
    return -1;

  for(int i=0;i<nc;i++)  {
    if (C[i]->ns > maxns)  {
      maxns = C[i]->ns;
      maxidx = i;
    }
  }
  return maxidx;

}

// deallocate the array of slices from allslices 
int dealloc_slice_array(slice **C, unsigned nc)  {
  for(int i=0;i<nc;i++)  {
    free(C[i]->S);
    free(C[i]);
  }
  free(C);
  return 0; 
}

// deep-copy slice
slice *copyslice(slice *slc)  {
  slice *newslc = (slice *)malloc(sizeof(slice));
  newslc->start = slc->start;
  newslc->stop = slc->stop;
  newslc->step = slc->step;
  newslc->ns = slc->ns;
  newslc->S = (unsigned *)malloc(newslc->ns*sizeof(unsigned));
  for(int i=0;i<slc->ns;i++)  
    newslc->S[i] = slc->S[i];
  return newslc;
}

// Algorithm 1: 
// repeat
//   take longest slice
//   modify Q
//   compute allslices
// 
// Note: must pass NULL for rptsc and slcs
int algorithm1(unsigned *X, int np, unsigned **rpoints, int *nrpt, slice ***slices, int *nslc)  {

  slice **C = NULL;
  unsigned nc = 0;

  // put slices in a slice list because we don't know
  // how many we'll end up with.
  slicenode *head = NULL;
  unsigned nn = 0;

  sequence *Q = initseq(X, np);
  while (1)  {

    int maxidx = -1;

    if (allslices(Q, &C, &nc, 0) != 0)  {
      break;
    }

    if ((maxidx = maxslice(C, nc)) == -1)  {
      break;
    }

    if (removeslice(Q, C[maxidx]) != 0)  {
      break;
    }

    // put slice in list
    insert_slice_list(copyslice(C[maxidx]), &head);
    nn++;

    // deallocate 
    dealloc_slice_array(C, nc);

  }

  // pack slices into an array
  slice **SS = (slice **)malloc(nn*sizeof(slice *));
  slicenode *s = head;
  for(int i=0;i<nn;i++)  {
    SS[i] = s->slice;
    s = s->nxt;
  }

  // deallocate the list
  delete_slice_list(head); 

  // delete everything in Q->D except Q->X, which is returned
  for(int i=0;i<Q->nd;i++)  {
    free(Q->D[i]);
  }
  free(Q->D);

  *rpoints = Q->X;     // points remaining
  *nrpt = Q->np;

  *slices = SS;     // slices found
  *nslc = nn;

  free(Q);
  return 0;
}



