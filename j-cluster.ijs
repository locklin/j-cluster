3 : 0''
if. UNAME-:'Linux' do.
  LIBCLUST=:  '/home/scott/src/j-cluster/libhcluster.so'
elseif. UNAME-:'Darwin' do.
  'platform not supported' 13!:8[10
elseif. UNAME-: 'Win' do.
  'platform not supported' 13!:8[10
elseif. do.
  'platform not supported' 13!:8[10
end.
)

NB. creates mask
isNum =: [: -. 128!:5 +. _ = ]

NB. Node* treecluster (int nrows, int ncolumns, double** data, int** mask,
NB.   double weight[], int transpose, char dist, char method, double** distmatrix);
treecluster =: 4 : 0
 'data wt dist meth distmx' =. y
 mask =. isNum data
 'nr nc' =. $ data
 cmd=. LIBCLUST,' treecluster * x x *d *x *d x c c *d'
 0 pick cmd cd nr;nc;data;mask;wt;0;dist;meth;distmx 
)


treetst =: 3 : 0
 dist=.'e' 
 meth =. 's'
 data =. y
 wt =. (}.$ y) $ 2.7 - 1.7
 mask =. 0 + isNum data
 'nr nc' =. $ data
 cmd=. LIBCLUST,' treeclusterj * x x *d *x *d x c c'
 0 pick cmd cd nr;nc;data;mask;wt;0;dist;meth 
)


treetst2 =: 3 : 0
 ('e';'s') treetst2 y
:
 'dist meth' =. x
 data =. y
 wt =. (}.$ y) $ 2.7 - 1.7
 mask =. 0 + isNum data
 'nr nc' =. $ data
 cmd=. LIBCLUST,' treeclusterj * x x *d *x *d x c c'
 0 pick cmd cd nr;nc;data;mask;wt;0;dist;meth 
)


mymetric =: 4 : 0
'i j'=.x
 mask =. 0 + isNum y
 n=. #y
 wt=. n$ 2.7 -2.7
 cmd=. LIBCLUST,' mymetric d c i *d *x *d i i'
0 pick cmd cd 'e';n;y;mask;wt;i;j
)


euclid =: 4 : 0
'i j'=.x
 mask =. 0 + isNum y
 n=. #y
 wt=. n$ 2.7 -2.7
 cmd=. LIBCLUST,' euclid d i *d *d *x *x *d i i i'
0 pick cmd cd n;y;y;mask;mask;wt;i;j;0
)


euclid2 =: 4 : 0
'i j'=.x
 mask =. 0 + isNum y
 n=. 4 NB.}.$ y
 wt=. n$ 2.7 -2.7
 cmd=. LIBCLUST,' euclid2 d *d i i i '
0 pick cmd cd (,y);n;i;j
)

cuttree =: 3 : 0
 'nelements tree nclust'=.y
 clustid =. nelements $ 0
 cmd=. LIBCLUST,' cuttree n x x x *i'
 4 pick cmd cd nelements;tree;nclust;clustid
)

clusterdistance =: 4 : 0
'data wt'
)

median=: 3 : 0
n=.#y
cmd=. LIBCLUST,' median d x *d'
0 pick cmd cd n;y
)

NB. use like 'e' distancematrix data;wts
distancematrix=: 4 : 0
'data wts' =. y
mask=. isNum data
'nr nc'=. $ data
cmd=. LIBCLUST,' distancematrix *d x x *d *x *d c x'
0 pick cmd cd nr;nc;data;mask;wts;x;0
)


somcluster =: 3 : 0
 cmd=. LIBCLUST,' somcluster n x x *d *x *d x x x d x c *d x'
)



NB. henry rich's implementation of  Levenshtein  distance from j-list
levdist=: 4 : 0"1
'a b'=. (/: #&>)x;y
z=. >: iz =. i.#b
for_j. a do.
   z=. <./\&.(-&iz) (>: <. (j ~: b) + |.!.j_index) z
end.
{:z
) 
NB. Henry Rich does Jaccard
tanimoto =: (+&#   %/@:-   2 1&*@:(+/@:e.))&~."1 

NB. JP Jacobs
euclidean=: +/&:*:@(-"1)/ 



NB. Further reduction in execution time can be gotten by using integer arithmetic:
NB. d=: +/&:*:@(-"1)/&([: <. (2^20) * ]) 

upTri=: , #~ [: , [: -.@>:/~@i. #

