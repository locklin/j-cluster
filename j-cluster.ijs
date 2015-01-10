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
masker =: 0+ [: -. 128!:5 +. _ = ]

unDoub =: 3 : 0
 (15!:14 <'y') +(8*{:$y)*i.{.$y
)

unDoub2 =: 3 : 0
 yy =. 2 ic ,y
 (15!:14 <'yy') +(4*{:$y)*i.{.$y
)


treetst2 =: 3 : 0
 ('e';'a') treetst2 y
:
 'dist meth' =. x
 wt =. ({: $ y) $ 0.5 NB.1.7 - 1.7
 mask =. ,masker2 y
 'nr nc' =. $ y
 cmd=. LIBCLUST,' treeclusterj * x x *x *x *d i c c'
 0 pick cmd cd nr;nc;(unDoub y);(unDoub2 mask);wt;0;dist;meth 
)

treetst =: 3 : 0
 ('e';'s') treetst y
:
 'dist meth' =. x
 wt =. ({: $ y) $ 2.7 - 1.7
 mask =.  masker y
 'nr nc' =. $ y
 cmd=. LIBCLUST,' treecluster * x x *x *x *d x c c *x'
 0 pick cmd cd nr;nc;(unDoub y); (unDoub2 mask) ;wt;0;dist;meth;''
)


mymetric =: 4 : 0
'i j'=.x
 mask =. masker y
 n=. #y
 wt=. n$ 3.7 -2.7
 cmd=. LIBCLUST,' mymetric d c i *x *x *d i i'
0 pick cmd cd 'e';n;(unDoub y);(unDoub mask);wt;i;j
)



euclid2 =: 4 : 0
'i j'=.x
 mask =. masker y
 n=. {:$ y 
 wt=. n$ 3.7 -2.7
 cmd=. LIBCLUST,' euclid2 d *x i i i '
0 pick cmd cd (unDoub y);n;i;j
)

cutree =: 3 : 0
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
 mask=. masker data
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

