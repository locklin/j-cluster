/******************************************************************************/
/* The C Clustering Library.
 * Copyright (C) 2002 Michiel Jan Laurens de Hoon.
 *
 * This library was written at the Laboratory of DNA Information Analysis,
 * Human Genome Center, Institute of Medical Science, University of Tokyo,
 * 4-6-1 Shirokanedai, Minato-ku, Tokyo 108-8639, Japan.
 * Contact: mdehoon 'AT' gsc.riken.jp
 * 
 * Permission to use, copy, modify, and distribute this software and its
 * documentation with or without modifications and for any purpose and
 * without fee is hereby granted, provided that any copyright notices
 * appear in all copies and that both those copyright notices and this
 * permission notice appear in supporting documentation, and that the
 * names of the contributors or copyright holders not be used in
 * advertising or publicity pertaining to distribution of the software
 * without specific prior permission.
 * 
 * THE CONTRIBUTORS AND COPYRIGHT HOLDERS OF THIS SOFTWARE DISCLAIM ALL
 * WARRANTIES WITH REGARD TO THIS SOFTWARE, INCLUDING ALL IMPLIED
 * WARRANTIES OF MERCHANTABILITY AND FITNESS, IN NO EVENT SHALL THE
 * CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY SPECIAL, INDIRECT
 * OR CONSEQUENTIAL DAMAGES OR ANY DAMAGES WHATSOEVER RESULTING FROM LOSS
 * OF USE, DATA OR PROFITS, WHETHER IN AN ACTION OF CONTRACT, NEGLIGENCE
 * OR OTHER TORTIOUS ACTION, ARISING OUT OF OR IN CONNECTION WITH THE USE
 * OR PERFORMANCE OF THIS SOFTWARE.
 * 
 */
#include <time.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>

#ifndef min
#define min(x, y)	((x) < (y) ? (x) : (y))
#endif
#ifndef max
#define	max(x, y)	((x) > (y) ? (x) : (y))
#endif

#ifdef WINDOWS
#  include <windows.h>
#endif

#define CLUSTERVERSION "2.0"

double** distancematrix (int ngenes, int ndata, char dist, double** data);
double farthest_distance(int n, double** dist);
double summed_distances(int n, double** dist);
int freedistmx(int nrow, double** distance);
void show_dists(int nrow,  double** distmx);

typedef struct {int left; int right; double distance;} Node;
/*
 * A Node struct describes a single node in a tree created by hierarchical
 * clustering. The tree can be represented by an array of n Node structs,
 * where n is the number of elements minus one. The integers left and right
 * in each Node struct refer to the two elements or subnodes that are joined
 * in this node. The original elements are numbered 0..nelements-1, and the
 * nodes -1..-(nelements-1). For each node, distance contains the distance
 * between the two subnodes that were joined.
 */
int freeNodes(Node* mynode);


Node* treecluster (char dist, char method, int nrows, int ncolumns, double** dst, double** data);

void cuttree (int nelements, Node* tree, int nclusters, int clusterid[]);

int dumpTree(int nc,Node* tree, double dist [], int lt[], int rt[]);


/* /\* Chapter 6 *\/ */
/* int pca(int m, int n, double** u, double** v, double* w); */

/* /\* Utility routines, currently undocumented *\/ */
 void sort(int n, const double data[], int index[]); 






