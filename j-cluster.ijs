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
 wt =. (}.$ y) # 2.7 - 1.7
 mask =. 0 + isNum data
 'nr nc' =. $ data
 cmd=. LIBCLUST,' treeclusterj * x x *d *x *d x c c'
 0 pick cmd cd nr;nc;data;mask;wt;0;dist;meth 
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

