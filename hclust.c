#include "hclust.h"
#include <stdio.h> 

/* ---------------------------------------------------------------------- */
int freedistmx(int n, double** dist)
{
  int i;
  for (i = 1; i < n; i++) {
    free(dist[i]);
  }
  free (dist);
  return 1;
}

/* ---------------------------------------------------------------------- */
int freeNodes(Node* mynodes)
{
  free(mynodes);
  return 1;
}
/* ---------------------------------------------------------------------- */



/* Print the data matrix */
void show_dists(int nrows, double** distMx){
  int i, j;
  printf("   Row: ");
  for(i=0; i<nrows-1; i++) printf("%9d", i);
  printf("\n");
  for(i=0; i<nrows; i++) {
    printf("Row %2d: ",i);
    for(j=0; j<i; j++) printf(" %f",distMx[i][j]);
    printf("\n");
  }
  printf("\n");
  return;
}


/*
static double euclid (int n, double** data1, double** data2, int index1, int index2,

Purpose
=======

The euclid routine calculates the weighted Euclidean distance between two
rows or columns in a matrix.

Arguments
=========

n      (input) int
The number of elements in a row or column. If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) double array
The data array containing the first vector.

data2  (input) double array
The data array containing the second vector.

weight (input) double[n]
The weights that are used to calculate the distance.

index1     (input) int
Index of the first row or column.

index2     (input) int
Index of the second row or column.


============================================================================
*/
static double euclid (int n, double** data1, double** data2, int index1, int index2) {
  double result = 0.;
  int i;
  for (i = 0; i < n; i++) {
	double term = data1[index1][i] - data2[index2][i];
        result += term*term;
    }
  return result;
}

/* ********************************************************************* */

/*

static double cityblock (int n, double** data1, double** data2, int index1, int index2)

Purpose
=======
The cityblock routine calculates the weighted "City Block" distance between
two rows or columns in a matrix. City Block distance is defined as the
absolute value of X1-X2 plus the absolute value of Y1-Y2 plus..., which is
equivalent to taking an "up and over" path.

Arguments (same as above)

============================================================================ */
static double cityblock (int n, double** data1, double** data2,int index1, int index2) {
  double result = 0.;
  int i;
  for (i = 0; i < n; i++) {
    double term = data1[index1][i] - data2[index2][i];
    result +=  fabs(term);
  }
  return result;
}

/* ********************************************************************* */

/*

static
double correlation (int n, double** data1, double** data2, int index1, int index2, 

Purpose
=======

The correlation routine calculates the weighted Pearson distance between two
rows or columns in a matrix. We define the Pearson distance as one minus the
Pearson correlation.
This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
(e.g., choose b = a + c).

Arguments
=========
(same as above)
============================================================================
*/
static double correlation (int n, double** data1, double** data2, int index1, 
			   int index2){
  double result = 0.;
  double sum1 = 0.;
  double sum2 = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  double tweight =0.;
  int i;
  for (i = 0; i < n; i++) {
    double term1 = data1[index1][i];
    double term2 = data2[index2][i];
    sum1 += term1;
    sum2 += term2;
    result += term1*term2;
    denom1 += term1*term1;
    denom2 += term2*term2;
    tweight += 1.0;
  }
  if(!tweight) return 0; 
  result -= sum1 * sum2 /tweight;
  denom1 -= sum1 * sum1/tweight;
  denom2 -= sum2 * sum2/tweight;
  if (denom1 <= 0) return 1; /* include '<' to deal with roundoff errors */
  if (denom2 <= 0) return 1; /* include '<' to deal with roundoff errors */
  result = result / sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}

/*
static double acorrelation (int n, double** data1, double** data2, int** mask1,
			    int** mask2, const double weight[], int index1, int index2, 
			    int transpose)
Purpose
=======

The acorrelation routine calculates the weighted Pearson distance between two
rows or columns, using the absolute value of the correlation.
This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
(e.g., choose b = a + c).

Arguments (same as above)

============================================================================
*/
static double acorrelation (int n, double** data1, double** data2, int index1, int index2) {
  double result = 0.;
  double sum1 = 0.;
  double sum2 = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  double tweight = 0.;
  int i;
  for (i = 0; i < n; i++) {
    double term1 = data1[index1][i];
    double term2 = data2[index2][i];
    sum1 += term1;
    sum2 += term2;
    result += term1*term2;
    denom1 += term1*term1;
    denom2 += term2*term2;
    tweight += 1.0;
  }
  result -= sum1 * sum2 / tweight;
  denom1 -= sum1 * sum1 / tweight;
  denom2 -= sum2 * sum2 / tweight;
  if (denom1 <= 0) return 1; /* include '<' to deal with roundoff errors */
  if (denom2 <= 0) return 1; /* include '<' to deal with roundoff errors */
  result = fabs(result) / sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}

/* ********************************************************************* */

/* static */
/* double ucorrelation (int n, double** data1, double** data2, int index1, int index2) */
/*
Purpose
=======

The ucorrelation routine calculates the weighted Pearson distance between two
rows or columns, using the uncentered version of the Pearson correlation. In the
uncentered Pearson correlation, a zero mean is used for both vectors even if
the actual mean is nonzero.
This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
(e.g., choose b = a + c).

Arguments (same as above)

============================================================================
*/
static double ucorrelation (int n, double** data1, double** data2, int index1, int index2){ 
  int i;
  double result = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  for (i = 0; i < n; i++)    { 
    double term1 = data1[index1][i];
    double term2 = data2[index2][i];
    result += term1*term2;
    denom1 += term1*term1;
    denom2 += term2*term2;
  }
  if (denom1==0.) return 1.;   /* these should be altered SCL */
  if (denom2==0.) return 1.;
  result = result / sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}

/* ********************************************************************* */

 /* static */
 /* double uacorrelation (int n, double** data1, double** data2, int index1, int index2)*/ 
 /*
   Purpose
=======

The uacorrelation routine calculates the weighted Pearson distance between two
rows or columns, using the absolute value of the uncentered version of the
Pearson correlation. In the uncentered Pearson correlation, a zero mean is used
for both vectors even if the actual mean is nonzero.
This definition yields a semi-metric: d(a,b) >= 0, and d(a,b) = 0 iff a = b.
but the triangular inequality d(a,b) + d(b,c) >= d(a,c) does not hold
(e.g., choose b = a + c).

Arguments (same as above)

============================================================================
*/
static double uacorrelation (int n, double** data1, double** data2, int index1, 
			     int index2) {
  int i;
  double result = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  
  for (i = 0; i < n; i++)    { 
    double term1 = data1[index1][i];
    double term2 = data2[index2][i];
    result += term1*term2;
    denom1 += term1*term1;
    denom2 += term2*term2;
  }
  
  if (denom1==0.) return 1.;
  if (denom2==0.) return 1.;
  result = fabs(result) / sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}

static const double* sortdata = NULL; /* used in the quicksort algorithm */

/* Helper function for sort. Previously, this was a nested function under
 * sort, which is not allowed under ANSI C.
 */
static int compare(const void* a, const void* b){ 
  const int i1 = *(const int*)a;
  const int i2 = *(const int*)b;
  const double term1 = sortdata[i1];
  const double term2 = sortdata[i2];
  if (term1 < term2) return -1;
  if (term1 > term2) return +1;
  return 0;
}

void sort(int n, const double data[], int index[])
/* Sets up an index table given the data, such that data[index[]] is in
 * increasing order. Sorting is done on the indices; the array data
 * is unchanged.
 */
{ int i;
  sortdata = data;
  for (i = 0; i < n; i++) index[i] = i;
  qsort(index, n, sizeof(int), compare);
}


static double* getrank (int n, double data[])
/* Calculates the ranks of the elements in the array data. Two elements with
 * the same value get the same rank, equal to the average of the ranks had the
 * elements different values. The ranks are returned as a newly allocated
 * array that should be freed by the calling routine. If getrank fails due to
 * a memory allocation error, it returns NULL.
 */
{ int i;
  double* rank;
  int* index;
  rank = malloc(n*sizeof(double));
  if (!rank) return NULL;
  index = malloc(n*sizeof(int));
  if (!index)
  { free(rank);
    return NULL;
  }
  /* Call sort to get an index table */
  sort (n, data, index);
  /* Build a rank table */
  for (i = 0; i < n; i++) rank[index[i]] = i;
  /* Fix for equal ranks */
  i = 0;
  while (i < n)
  { int m;
    double value = data[index[i]];
    int j = i + 1;
    while (j < n && data[index[j]] == value) j++;
    m = j - i; /* number of equal ranks found */
    value = rank[index[i]] + (m-1)/2.;
    for (j = i; j < i + m; j++) rank[index[j]] = value;
    i += m;
  }
  free (index);
  return rank;
}

/* *********************************************************************  */

/* static */
/* double spearman (int n, double** data1, double** data2, */
/*    int index1, int index2) */
/*
Purpose
=======

The spearman routine calculates the Spearman distance between two rows or
columns. The Spearman distance is defined as one minus the Spearman rank
correlation.

Arguments
=========

n      (input) int
The number of elements in a row . If transpose==0, then n is the number
of columns; otherwise, n is the number of rows.

data1  (input) double array
The data array containing the first vector.

data2  (input) double array
The data array containing the second vector.

index1     (input) int
Index of the first row .

index2     (input) int
Index of the second row .

============================================================================
*/
static double spearman (int n, double** data1, double** data2, int index1, 
			int index2){
  int i;
  int m = 0;
  double* rank1;
  double* rank2;
  double result = 0.;
  double denom1 = 0.;
  double denom2 = 0.;
  double avgrank;
  double* tdata1;
  double* tdata2;
  tdata1 = malloc(n*sizeof(double));
  if(!tdata1) return 0.0;                        /* Memory allocation error */
  tdata2 = malloc(n*sizeof(double));
  if(!tdata2) {                                     /* Memory allocation error */  
    free(tdata1);
    return 0.0;
  }
  for (i = 0; i < n; i++) {
      tdata1[m] = data1[index1][i];
      tdata2[m] = data2[index2][i];
      m++;
  }
  if (m==0) {
    free(tdata1);
    free(tdata2);
    return 0;
  }
  rank1 = getrank(m, tdata1);
  free(tdata1);
  if(!rank1) { 
    free(tdata2);
    return 0.0; /* Memory allocation error */
  }
  rank2 = getrank(m, tdata2);
  free(tdata2);
  if(!rank2) { /* Memory allocation error */
   free(rank1);
    return 0.0;
  }
  avgrank = 0.5*(m-1); /* Average rank */
  for (i = 0; i < m; i++) {
    const double value1 = rank1[i];
    const double value2 = rank2[i];
    result += value1 * value2;
    denom1 += value1 * value1;
    denom2 += value2 * value2;
  }
  /* Note: denom1 and denom2 cannot be calculated directly from the number
   * of elements. If two elements have the same rank, the squared sum of
   * their ranks will change.
   */
  free(rank1);
  free(rank2);
  result /= m;
  denom1 /= m;
  denom2 /= m;
  result -= avgrank * avgrank;
  denom1 -= avgrank * avgrank;
  denom2 -= avgrank * avgrank;
  if (denom1 <= 0) return 1; /* include '<' to deal with roundoff errors */
  if (denom2 <= 0) return 1; /* include '<' to deal with roundoff errors */
  result = result / sqrt(denom1*denom2);
  result = 1. - result;
  return result;
}


/* *********************************************************************  */

/*  static double kendall (int n, double** data1, double** data2, int index1, */
/*    int index2) */
/* Purpose
=======

The kendall routine calculates the Kendall distance between two
rows s. The Kendall distance is defined as one minus Kendall's tau.

Arguments (same as elsewhere)

============================================================================
*/
static double kendall (int n, double** data1, double** data2, int index1, int index2){
  int con = 0;
  int dis = 0;
  int exx = 0;
  int exy = 0;
  double denomx;
  double denomy;
  double tau;
  int i, j;
  for (i = 0; i < n; i++) {
    for (j = 0; j < i; j++) {
      double x1 = data1[index1][i];
      double x2 = data1[index1][j];
      double y1 = data2[index2][i];
      double y2 = data2[index2][j];
      if (x1 < x2 && y1 < y2) con++;
      if (x1 > x2 && y1 > y2) con++;
      if (x1 < x2 && y1 > y2) dis++;
      if (x1 > x2 && y1 < y2) dis++;
      if (x1 == x2 && y1 != y2) exx++;
      if (x1 != x2 && y1 == y2) exy++;
    }
  }
  denomx = con + dis + exx;
  denomy = con + dis + exy;
  if (denomx==0) return 1;
  if (denomy==0) return 1;
  tau = (con-dis)/sqrt(denomx*denomy);
  return 1.-tau;
}

/* *********************************************************************  */
/*
static double angle(int n, double** data1, double** data2,int index1, int index2)

Purpose
=======

The angle routine calculates the angular distance between two
rows s; aka angle(x,y) = arccos(x . y /(|x| |y|) )

Arguments (same as elsewhere)

*/

static double angle(int n, double** data1, double** data2, int index1, int index2){
  double ang;
  double xx = 0;
  double yy = 0;
  double xdoty = 0;
  int i;
  for (i = 0; i < n; i++) {
    xx +=  data1[index1][i] * data1[index1][i] ;
    yy +=  data2[index2][i] * data2[index2][i] ;
    xdoty += data1[index1][i] * data2[index2][i];
  }
  ang = acos( (xdoty/(xx*yy)));
  return ang;
}


/* *********************************************************************  */
/*
static double cosine(int n, double** data1, double** data2, int** mask1, int** mask2,
		    const double weight[],int index1, int index2, int transpose){
Purpose
=======

The cosine routine calculates the vector cosine between two
rows s; aka cosine(x,y) = x . y /(|x| |y|) 

Arguments (same as elsewhere)
*/

static double cosine(int n, double** data1, double** data2, int index1, int index2){
  double ang;
  double xx = 0;
  double yy = 0;
  double xdoty = 0;
  int i;
  for (i = 0; i < n; i++) {
    xx +=  data1[index1][i] * data1[index1][i] ;
    yy +=  data2[index2][i] * data2[index2][i] ;
    xdoty += data1[index1][i] * data2[index2][i];
  }
  ang =  xdoty / (xx*yy);
  return ang;
}


/* *********************************************************************  */
/*
static double chebyshev(int n, double** data1, double** data2,int index1,int index2) 

Purpose
=======

The chebyshev distance calculates the vector cosine between two
rows; aka cosine(x,y) = x . y /(|x| |y|) 

Arguments (same as elsewhere)
*/

static double chebyshev(int n, double** data1, double** data2, 
			int index1, int index2){
  double linf = 0;
  int i;
  for (i = 0; i < n; i++) {
      linf = max(linf, abs( data1[index1][i] - data2[index2][i]));
  }
  return linf;
}





/* *********************************************************************  */

static double(*setmetric(char dist))  (int, double**, double**, int, int){ 
  switch(dist)  { 
  case 'e': return &euclid;
  case 'b': return &cityblock;
  case 'c': return &correlation;
  case 'a': return &acorrelation;
  case 'u': return &ucorrelation;
  case 'x': return &uacorrelation;
  case 's': return &spearman;
  case 'k': return &kendall;
  case 'o': return &angle; 
  case 'n': return &cosine; 
  case 'y': return &chebyshev;
  default: return &euclid;
  }
  return NULL; /* Never get here */
}



/*
double** distancematrix (int nrows, int ncolumns, double** data,
  int** mask, double weights[], char dist, int transpose)

Purpose
=======

The distancematrix routine calculates the distance matrix between genes or
microarrays using their measured gene expression data. Several distance measures
can be used. The routine returns a pointer to a ragged array containing the
distances between the genes. As the distance matrix is symmetric, with zeros on
the diagonal, only the lower triangular half of the distance matrix is saved.
The distancematrix routine allocates space for the distance matrix. If the
parameter transpose is set to a nonzero value, the distances between the columns
(microarrays) are calculated, otherwise distances between the rows (genes) are
calculated.
If sufficient space in memory cannot be allocated to store the distance matrix,
the routine returns a NULL pointer, and all memory allocated so far for the
distance matrix is freed.


Arguments
=========

nrows      (input) int
The number of rows in the gene expression data matrix (i.e., the number of
genes)

ncolumns   (input) int
The number of columns in the gene expression data matrix (i.e., the number of
microarrays)

data       (input) double[nrows][ncolumns]
The array containing the gene expression data.

dist       (input) char
Defines which distance measure is used, as given by the table:
dist=='e': Euclidean distance
dist=='b': City-block distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.


========================================================================
*/
double** distancematrix (int nrows, int ncolumns, char dist, double** data){ 
  /* First determine the size of the distance matrix */
  int i,j;
  double** matrix;
  /* Set the metric function as indicated by dist */
  double (*metric)  (int, double**, double**, int, int) =  setmetric(dist);
  if (nrows < 2) return NULL;
  /* Set up the ragged array */
  matrix = malloc(nrows*sizeof(double*));
  if(matrix==NULL) return NULL; /* Not enough memory available */
  matrix[0] = NULL;
  /* The zeroth row has zero columns. We allocate it anyway for convenience.*/
  for (i = 1; i < nrows; i++){ 
    matrix[i] = malloc(i*sizeof(double));
    if (matrix[i]==NULL) break; /* Not enough memory available */
  }
  if (i < nrows){ /* break condition encountered */
    j = i;
    for (i = 1; i < j; i++) free(matrix[i]);
    return NULL;
  }

  for (i = 1; i < nrows; i++) {
    for (j = 0; j < i; j++) {
       matrix[i][j]=metric(ncolumns,data,data,i,j);
    }
  }
  return matrix;
}

/*
static double find_closest_pair(int n, double** distmatrix, int* ip, int* jp)

This function searches the distance matrix to find the pair with the shortest
distance between them. The indices of the pair are returned in ip and jp; the
distance itself is returned by the function.

n          (input) int
The number of elements in the distance matrix.

distmatrix (input) double**
A ragged array containing the distance matrix. The number of columns in each
row is one less than the row index.

ip         (output) int*
A pointer to the integer that is to receive the first index of the pair with
the shortest distance.

jp         (output) int*
A pointer to the integer that is to receive the second index of the pair with
the shortest distance.
*/
static double find_closest_pair(int n, double** distmatrix, int* ip, int* jp) { 
  int i, j;
  double temp;
  double distance = distmatrix[1][0];
  *ip = 1;
  *jp = 0;
  for (i = 1; i < n; i++){ 
    for (j = 0; j < i; j++){ 
      temp = distmatrix[i][j];
      if (temp<distance) { 
	distance = temp;
        *ip = i;
        *jp = j;
      }
    }
  }
  return distance;
}

/* these two functions are useful for geometers and topologists */
double farthest_distance(int n, double** distmatrix) {
  int i, j;
  double temp;
  double distance = distmatrix[1][0];
  for (i = 1; i < n; i++){ 
    for (j = 0; j < i; j++){ 
      temp = distmatrix[i][j];
      if (temp>distance) { 
	distance = temp;
      }
    }
  }
  return distance;
}


 double summed_distances(int n, double** distmatrix) {
  int i, j;
  double distance = distmatrix[1][0];
  for (i = 1; i < n; i++){ 
    for (j = 0; j < i; j++){ 
      distance = distance + distmatrix[i][j];
    }
  }
  return distance;
}


/*
void cuttree (int nelements, Node* tree, int nclusters, int clusterid[])  

Purpose
=======

The cuttree routine takes the output of a hierarchical clustering routine, and
divides the elements in the tree structure into clusters based on the
hierarchical clustering result. The number of clusters is specified by the user.

Arguments
=========

nelements      (input) int
The number of elements that were clustered.

tree           (input) Node[nelements-1]
The clustering solution. Each node in the array describes one linking event,
with tree[i].left and tree[i].right representig the elements that were joined.
The original elements are numbered 0..nelements-1, nodes are numbered
-1..-(nelements-1).

nclusters      (input) int
The number of clusters to be formed.

clusterid      (output) int[nelements]
The number of the cluster to which each element was assigned. Space for this
array should be allocated before calling the cuttree routine. If a memory
error occured, all elements in clusterid are set to -1.

========================================================================
*/
void cuttree (int nelements, Node* tree, int nclusters, int clusterid[])  { 
  int i, j, k;
  int icluster = 0;
  const int n = nelements-nclusters; /* number of nodes to join */
  int* nodeid;
  for (i = nelements-2; i >= n; i--)
  { k = tree[i].left;
    if (k>=0)
    { clusterid[k] = icluster;
      icluster++;
    }
    k = tree[i].right;
    if (k>=0)
    { clusterid[k] = icluster;
      icluster++;
    }
  }
  nodeid = malloc(n*sizeof(int));
  if(!nodeid)
  { for (i = 0; i < nelements; i++) clusterid[i] = -1;
    return;
  }
  for (i = 0; i < n; i++) nodeid[i] = -1;
  for (i = n-1; i >= 0; i--)
  { if(nodeid[i]<0) 
    { j = icluster;
      nodeid[i] = j;
      icluster++;
    }
    else j = nodeid[i];
    k = tree[i].left;
    if (k<0) nodeid[-k-1] = j; else clusterid[k] = j;
    k = tree[i].right;
    if (k<0) nodeid[-k-1] = j; else clusterid[k] = j;
  }
  free(nodeid);
  return;
}

/* ******************************************************************** */

int dumpTree(int nr, Node* tree, double dist [], int lt[], int rt[]) {
  int i;
  for (i = 0; i < nr; i++) {
    lt[i] = tree[i].left;
    rt[i] = tree[i].right;
    dist[i] = tree[i].distance;
  }
  return 1;
} 


/* ******************************************************************** */

/*
static
Node* pclcluster (int nrows, int ncolumns, double** data, double** distmatrix, 
char dist, int transpose)

Purpose
=======

The pclcluster routine performs clustering using pairwise centroid-linking
on a given set of gene expression data, using the distance metric given by dist.

Arguments
=========

nrows     (input) int
The number of rows in the gene expression data matrix, equal to the number of
genes.

ncolumns  (input) int
The number of columns in the gene expression data matrix, equal to the number of
microarrays.

data       (input) double[nrows][ncolumns]
The array containing the gene expression data.

dist       (input) char
Defines which distance measure is used, as given by the table:
dist=='e': Euclidean distance
dist=='b': City-block distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.

distmatrix (input) double**
The distance matrix. This matrix is precalculated by the calling routine
treecluster. The pclcluster routine modifies the contents of distmatrix, but
does not deallocate it.

Return value
============

A pointer to a newly allocated array of Node structs, describing the
hierarchical clustering solution consisting of nelements-1 nodes. Depending on
whether genes (rows) or microarrays (columns) were clustered, nelements is
equal to nrows or ncolumns. See src/cluster.h for a description of the Node
structure.
If a memory error occurs, pclcluster returns NULL.
========================================================================
*/
static Node* pclcluster (int nrows, int ncolumns, double** data, 
			 double** distmatrix, char dist){ 
  int i, j;
  const int nelements = nrows;
  int inode;
  const int ndata =  ncolumns;
  const int nnodes = nelements - 1;

  /* Set the metric function as indicated by dist */
  double (*metric)  (int, double**, double**, int, int) = setmetric(dist);
  Node* result;
  double** newdata;
  int* distid = malloc(nelements*sizeof(int));
  if(!distid) return NULL;
  result = malloc(nnodes*sizeof(Node));
  if(!result)  { 
    free(distid);
    return NULL;
  }

  newdata = malloc(nrows*sizeof(double*));
  for(i = 0; i< nrows; i++) {
    newdata[i] = malloc(ncolumns*sizeof(double));
  }

  for (i = 0; i < nelements; i++) distid[i] = i;
  /* To remember which row/column in the distance matrix contains what */

  /* Storage for node data */
  for (i = 0; i < nelements; i++){
    for (j = 0; j < ndata; j++){
      newdata[i][j] = data[j][i];
    }
  }
  data = newdata;
  

  for (inode = 0; inode < nnodes; inode++){ 
    /* Find the pair with the shortest distance */
    int is = 1;
    int js = 0;
    result[inode].distance = find_closest_pair(nelements-inode, distmatrix, &is, &js);
    result[inode].left = distid[js];
    result[inode].right = distid[is];
    /* Make node js the new node */
    for (i = 0; i < ndata; i++){
      data[js][i] = data[js][i] + data[is][i];
    }
    free(data[is]);
    data[is] = data[nnodes-inode];
  
    /* Fix the distances */
    distid[is] = distid[nnodes-inode];
    for (i = 0; i < is; i++)
      distmatrix[is][i] = distmatrix[nnodes-inode][i];
    for (i = is + 1; i < nnodes-inode; i++)
      distmatrix[i][is] = distmatrix[nnodes-inode][i];

    distid[js] = -inode-1;
    for (i = 0; i < js; i++)
      distmatrix[js][i] = metric(ndata,data,data,js,i);
    for (i = js + 1; i < nnodes-inode; i++)
      distmatrix[i][js] = metric(ndata,data,data,js,i);
  }

  /* Free temporarily allocated space */
  free(data[0]);
  free(data);
  free(distid);
  return result;
}

/* ******************************************************************** */

static
int nodecompare(const void* a, const void* b)
/* Helper function for qsort. */
{ const Node* node1 = (const Node*)a;
  const Node* node2 = (const Node*)b;
  const double term1 = node1->distance;
  const double term2 = node2->distance;
  if (term1 < term2) return -1;
  if (term1 > term2) return +1;
  return 0;
}

/* ---------------------------------------------------------------------- */


/*
static
Node* pslcluster (int nrows, int ncolumns, double** data, int** mask,
  double weight[], double** distmatrix, char dist, int transpose)

Purpose
=======

The pslcluster routine performs single-linkage hierarchical clustering, using
either the distance matrix directly, if available, or by calculating the
distances from the data array. This implementation is based on the SLINK
algorithm, described in:
Sibson, R. (1973). SLINK: An optimally efficient algorithm for the single-link
cluster method. The Computer Journal, 16(1): 30-34.
The output of this algorithm is identical to conventional single-linkage
hierarchical clustering, but is much more memory-efficient and faster. Hence,
it can be applied to large data sets, for which the conventional single-
linkage algorithm fails due to lack of memory.


Arguments
=========

nrows     (input) int
The number of rows in the gene expression data matrix, equal to the number of
genes.

ncolumns  (input) int
The number of columns in the gene expression data matrix, equal to the number of
microarrays.

data       (input) double[nrows][ncolumns]
The array containing the gene expression data.

mask       (input) int[nrows][ncolumns]
This array shows which data values are missing. If
mask[i][j] == 0, then data[i][j] is missing.

weight (input) double[n]
The weights that are used to calculate the distance. The length of this vector
is ncolumns if genes are being clustered, and nrows if microarrays are being
clustered.

transpose  (input) int
If transpose==0, the rows of the matrix are clustered. Otherwise, columns
of the matrix are clustered.

dist       (input) char
Defines which distance measure is used, as given by the table:
dist=='e': Euclidean distance
dist=='b': City-block distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.

distmatrix (input) double**
The distance matrix. If the distance matrix is passed by the calling routine
treecluster, it is used by pslcluster to speed up the clustering calculation.
The pslcluster routine does not modify the contents of distmatrix, and does
not deallocate it. If distmatrix is NULL, the pairwise distances are calculated
by the pslcluster routine from the gene expression data (the data and mask
arrays) and stored in temporary arrays. If distmatrix is passed, the original
gene expression data (specified by the data and mask arguments) are not needed
and are therefore ignored.


Return value
============

A pointer to a newly allocated array of Node structs, describing the
hierarchical clustering solution consisting of nelements-1 nodes. Depending on
whether genes (rows) or microarrays (columns) were clustered, nelements is
equal to nrows or ncolumns. See src/cluster.h for a description of the Node
structure.
If a memory error occurs, pslcluster returns NULL.

========================================================================
*/
static Node* pslcluster (int nrows, int ncolumns, double** data, 
			 double** distmatrix, char dist) { 
  int i, j, k;
  const int nelements = nrows; 
  const int nnodes = nelements - 1;
  int* vector;
  double* temp;
  int* index;
  Node* result;
  temp = malloc(nnodes*sizeof(double));
  if(!temp) return NULL;
  index = malloc(nelements*sizeof(int));
  if(!index)
    { free(temp);
      return NULL;
    }
  vector = malloc(nnodes*sizeof(int));
  if(!vector)
    { free(index);
      free(temp);
      return NULL;
    }
  result = malloc(nelements*sizeof(Node));
  if(!result){
    free(vector);
    free(index);
    free(temp);
    return NULL;
  }

  for (i = 0; i < nnodes; i++) vector[i] = i;

  if(distmatrix){ 
    for (i = 0; i < nrows; i++){ 
      result[i].distance = DBL_MAX;
      for (j = 0; j < i; j++) temp[j] = distmatrix[i][j];
      for (j = 0; j < i; j++){
       k = vector[j];
       if (result[j].distance >= temp[j]){
         if (result[j].distance < temp[k]) temp[k] = result[j].distance;
	 result[j].distance = temp[j];
	 vector[j] = i;
       } else if (temp[j] < temp[k]) temp[k] = temp[j];
      }
      for (j = 0; j < i; j++){
        if (result[j].distance >= result[vector[j]].distance) vector[j] = i;
      }
    }
  } else  { 
    const int ndata = ncolumns; /*nrows;*/
    /* Set the metric function as indicated by dist */
    double (*metric)  (int, double**, double**, int, int) =  setmetric(dist);
    for (i = 0; i < nelements; i++){
      result[i].distance = DBL_MAX;
      for (j = 0; j < i; j++) temp[j] =
        metric(ndata, data, data,  i, j);
      for (j = 0; j < i; j++){
	k = vector[j];
        if (result[j].distance >= temp[j]){
	  if (result[j].distance < temp[k]) temp[k] = result[j].distance;
          result[j].distance = temp[j];
          vector[j] = i;
        }
        else if (temp[j] < temp[k]) temp[k] = temp[j];
      }
      for (j = 0; j < i; j++)
        if (result[j].distance >= result[vector[j]].distance) vector[j] = i;
    }
  }
  free(temp);

  for (i = 0; i < nnodes; i++) result[i].left = i;
  qsort(result, nnodes, sizeof(Node), nodecompare);

  for (i = 0; i < nelements; i++) index[i] = i;
  for (i = 0; i < nnodes; i++){
    j = result[i].left;
    k = vector[j];
    result[i].left = index[j];
    result[i].right = index[k];
    index[k] = -i-1;
  }
  free(vector);
  free(index);

  result = realloc(result, nnodes*sizeof(Node));

  return result;
}
/* ******************************************************************** */


/*
static Node* pmlcluster (int nelements, double** distmatrix)

Purpose
=======

The pmlcluster routine performs clustering using pairwise maximum- (complete-)
linking on the given distance matrix.

Arguments
=========

nelements     (input) int
The number of elements to be clustered.

distmatrix (input) double**
The distance matrix, with nelements rows, each row being filled up to the
diagonal. The elements on the diagonal are not used, as they are assumed to be
zero. The distance matrix will be modified by this routine.

Return value
============

A pointer to a newly allocated array of Node structs, describing the
hierarchical clustering solution consisting of nelements-1 nodes. Depending on
whether genes (rows) or microarrays (columns) were clustered, nelements is
equal to nrows or ncolumns. See src/cluster.h for a description of the Node
structure.
If a memory error occurs, pmlcluster returns NULL.
========================================================================
*/
static Node* pmlcluster (int nelements, double** distmatrix) { 
  int j,n;
  int* clusterid;
  Node* result;

  clusterid = malloc(nelements*sizeof(int));
  if(!clusterid) return NULL;
  result = malloc((nelements-1)*sizeof(Node));
  if (!result)
  { free(clusterid);
    return NULL;
  }

  /* Setup a list specifying to which cluster a gene belongs */
  for (j = 0; j < nelements; j++) clusterid[j] = j;

  for (n = nelements; n > 1; n--)
  { int is = 1;
    int js = 0;
    result[nelements-n].distance = find_closest_pair(n, distmatrix, &is, &js);

    /* Fix the distances */
    for (j = 0; j < js; j++)
      distmatrix[js][j] = max(distmatrix[is][j],distmatrix[js][j]);
    for (j = js+1; j < is; j++)
      distmatrix[j][js] = max(distmatrix[is][j],distmatrix[j][js]);
    for (j = is+1; j < n; j++)
      distmatrix[j][js] = max(distmatrix[j][is],distmatrix[j][js]);

    for (j = 0; j < is; j++) distmatrix[is][j] = distmatrix[n-1][j];
    for (j = is+1; j < n-1; j++) distmatrix[j][is] = distmatrix[n-1][j];

    /* Update clusterids */
    result[nelements-n].left = clusterid[is];
    result[nelements-n].right = clusterid[js];
    clusterid[js] = n-nelements-1;
    clusterid[is] = clusterid[n-1];
  }
  free(clusterid);

  return result;
}

/* ******************************************************************* */


/*
static Node* palcluster (int nelements, double** distmatrix)

Purpose
=======

The palcluster routine performs clustering using pairwise average
linking on the given distance matrix.

Arguments
=========

nelements     (input) int
The number of elements to be clustered.

distmatrix (input) double**
The distance matrix, with nelements rows, each row being filled up to the
diagonal. The elements on the diagonal are not used, as they are assumed to be
zero. The distance matrix will be modified by this routine.

Return value
============

A pointer to a newly allocated array of Node structs, describing the
hierarchical clustering solution consisting of nelements-1 nodes. Depending on
whether genes (rows) or microarrays (columns) were clustered, nelements is
equal to nrows or ncolumns. See src/cluster.h for a description of the Node
structure.
If a memory error occurs, palcluster returns NULL.
========================================================================
*/
static Node* palcluster (int nelements, double** distmatrix){ 
  int j,n;
  int* clusterid;
  int* number;
  Node* result;

  clusterid = malloc(nelements*sizeof(int));
  if(!clusterid) return NULL;
  number = malloc(nelements*sizeof(int));
  if(!number)  { 
    free(clusterid);
    return NULL;
  }
  result = malloc((nelements-1)*sizeof(Node));
  if (!result)  { 
    free(clusterid);
    free(number);
    return NULL;
  }

  /* Setup a list specifying to which cluster a gene belongs, and keep track
   * of the number of elements in each cluster (needed to calculate the
   * average). */
  for (j = 0; j < nelements; j++) {
    number[j] = 1;
    clusterid[j] = j;
  }

  for (n = nelements; n > 1; n--) {
    int sum;
    int is = 1;
    int js = 0;
    result[nelements-n].distance = find_closest_pair(n, distmatrix, &is, &js);

    /* Save result */
    result[nelements-n].left = clusterid[is];
    result[nelements-n].right = clusterid[js];

    /* Fix the distances */
    sum = number[is] + number[js];
    for (j = 0; j < js; j++) {
      distmatrix[js][j] = distmatrix[is][j]*number[is]
	+ distmatrix[js][j]*number[js];
      distmatrix[js][j] /= sum;
    }
    for (j = js+1; j < is; j++) {
      distmatrix[j][js] = distmatrix[is][j]*number[is]
                        + distmatrix[j][js]*number[js];
      distmatrix[j][js] /= sum;
    }
    for (j = is+1; j < n; j++) {
      distmatrix[j][js] = distmatrix[j][is]*number[is]
                        + distmatrix[j][js]*number[js];
      distmatrix[j][js] /= sum;
    }

    for (j = 0; j < is; j++) distmatrix[is][j] = distmatrix[n-1][j];
    for (j = is+1; j < n-1; j++) distmatrix[j][is] = distmatrix[n-1][j];

    /* Update number of elements in the clusters */
    number[js] = sum;
    number[is] = number[n-1];

    /* Update clusterids */
    clusterid[js] = n-nelements-1;
    clusterid[is] = clusterid[n-1];
  }
  free(clusterid);
  free(number);

  return result;
}

/* ******************************************************************* */


/* 
Node* treecluster (int nrows, int ncolumns, double** data, int** mask,
		   double weight[], int transpose, char dist, char method, 
		   double** distmatrix) 


Purpose
=======

The treecluster routine performs hierarchical clustering using pairwise
single-, maximum-, centroid-, or average-linkage, as defined by method, on a
given set of gene expression data, using the distance metric given by dist.
If successful, the function returns a pointer to a newly allocated Tree struct
containing the hierarchical clustering solution, and NULL if a memory error
occurs. The pointer should be freed by the calling routine to prevent memory
leaks.

Arguments
=========

nrows     (input) int
The number of rows in the data matrix, equal to the number of genes.

ncolumns  (input) int
The number of columns in the data matrix, equal to the number of microarrays.

data       (input) double[nrows][ncolumns]
The array containing the data of the vectors to be clustered.

mask       (input) int[nrows][ncolumns]
This array shows which data values are missing. If mask[i][j]==0, then
data[i][j] is missing.

weight (input) double array[n]
The weights that are used to calculate the distance.

transpose  (input) int
If transpose==0, the rows of the matrix are clustered. Otherwise, columns
of the matrix are clustered.

dist       (input) char
Defines which distance measure is used, as given by the table:
dist=='e': Euclidean distance
dist=='b': City-block distance
dist=='c': correlation
dist=='a': absolute value of the correlation
dist=='u': uncentered correlation
dist=='x': absolute uncentered correlation
dist=='s': Spearman's rank correlation
dist=='k': Kendall's tau
For other values of dist, the default (Euclidean distance) is used.

method     (input) char
Defines which hierarchical clustering method is used:
method=='s': pairwise single-linkage clustering
method=='m': pairwise maximum- (or complete-) linkage clustering
method=='a': pairwise average-linkage clustering
method=='c': pairwise centroid-linkage clustering
For the first three, either the distance matrix or the gene expression data is
sufficient to perform the clustering algorithm. For pairwise centroid-linkage
clustering, however, the gene expression data are always needed, even if the
distance matrix itself is available.

distmatrix (input) double**
The distance matrix. If the distance matrix is zero initially, the distance
matrix will be allocated and calculated from the data by treecluster, and
deallocated before treecluster returns. If the distance matrix is passed by the
calling routine, treecluster will modify the contents of the distance matrix as
part of the clustering algorithm, but will not deallocate it. The calling
routine should deallocate the distance matrix after the return from treecluster.

Return value
============

A pointer to a newly allocated array of Node structs, describing the
hierarchical clustering solution consisting of nelements-1 nodes. Depending on
whether genes (rows) or microarrays (columns) were clustered, nelements is
equal to nrows or ncolumns. See src/cluster.h for a description of the Node
structure.
If a memory error occurs, treecluster returns NULL.

========================================================================
*/
Node* treecluster (char dist, char method, int nrows, int ncolumns,double** distmatrix, double** data) {
  Node* result = NULL;
  const int nelements = nrows;
  const int ldistmatrix = (distmatrix==NULL && method!='s') ? 1 : 0;
  if (nelements < 2) return NULL;
  /* Calculate the distance matrix if the user didn't give it */
  if(ldistmatrix)  { 
    distmatrix =
      distancematrix(nrows, ncolumns, dist, data);
    if (!distmatrix) return NULL; /* Insufficient memory */
  }

  switch(method) {
  case 's':
    result = pslcluster(nrows, ncolumns, data, distmatrix, dist);
    break;
  case 'm':
    result = pmlcluster(nelements, distmatrix);
    break;
  case 'a':
    result = palcluster(nelements, distmatrix);
    break;
  case 'c':
    result = pclcluster(nrows, ncolumns, data, distmatrix, dist);
    break;
  }

  if(ldistmatrix) { 
    freedistmx(nelements,distmatrix);
  }
  return result;
}
