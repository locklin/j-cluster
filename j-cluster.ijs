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

NB. available distance metrics
NB. case 'e': return &euclid;
NB. case 'b': return &cityblock;
NB. case 'c': return &correlation;
NB. case 'a': return &acorrelation;
NB. case 'u': return &ucorrelation;
NB. case 'x': return &uacorrelation;
NB. case 's': return &spearman;
NB. case 'k': return &kendall;
NB. case 'o': return &angle; 
NB. case 'n': return &cosine; 
NB. case 'y': return &chebyshev;

NB. creates mask
masker =: 0+ [: -. 128!:5 +. _ = ]

unDoub =: 3 : 0
 (15!:14 <'y') +(8*{:$y)*i.{.$y
)

unInt =: 3 : 0
 yy =. 2 ic ,y
 (15!:14 <'yy') +(4*{:$y)*i.{.$y
)


treetst =: 3 : 0
 ('e';'s') treetst y
:
 'dist meth' =. x
 wt =. ({: $ y) $ 2.7 - 1.7
 mask =.  masker y
 'nr nc' =. $ y
 cmd=. LIBCLUST,' treecluster * x x *x *x *d x c c *x'
 0 pick cmd cd nr;nc;(unDoub y); (unInt mask) ;wt;0;dist;meth;(<<0)
)


cutree =: 3 : 0
 'nelements tree nclust'=.y
 clustid =. nelements $ 0
 cmd=. LIBCLUST,' cuttree n x x x *i'
 4 pick cmd cd nelements;tree;nclust;clustid
)

median=: 3 : 0
n=.#y
cmd=. LIBCLUST,' median d x *d'
0 pick cmd cd n;y
)


NB. this should probably never be used.
distancematrix=: 3 : 0
 'e' distancematrix y
:
 wts =. ({: $ y) $ 2.7 - 1.7
 mask=. masker y
 'nr nc' =. $ y
 mask =. masker y
 cmd=. LIBCLUST,' distancematrix * x x *x *x *d c x'
 0 pick cmd cd nr;nc;(unDoub y);(unInt mask);wts;x;0
)


NB. probably not needed
clusterdistance =: 4 : 0
 'nr nc' =. $ y
 wts =. ({: $ y) $ 2.7 - 1.7
 mask=. masker y
 cmd =. LIBCLUST,' clusterdistance d i i *x *x *d i i *i *i c c i'
)

NB. probably also not needed
somcluster =: 3 : 0
 cmd=. LIBCLUST,' somcluster n x x *d *x *d x x x d x c *d x'
)


NB. useful for comparison to lapack to see what your overhead is like
pcaC =: 3 : 0
 'nr nc' =. $ y
 nn =. (nr <. nc) 
 u =. y
 v =. (nn,nn) $ 1.1 - 1.1
 w =. (nn,nn) $ 1.1 - 1.1
 wts =. ({: $ y) $ 2.7 - 1.7
 cmd =. LIBCLUST,' pca i i i *x *x *x'
 cmd cd nr;nc;(unDoub u);(unDoub v);(unDoub w)
 u;v;w
) 


calcweights=: 4 : 0
 'kind cut exp' =. x
 'nr nc' =. $ y
 wts =. ({: $ y) $ 2.7 - 1.7
 mask=. masker y
 res =. nr $ 2.2 - 2.2
 cmd =. LIBCLUST,' calculate_weights *d i i *x *x *d *d i c d d'
 6 pick cmd cd nr;nc;(unDoub y);(unInt mask);wts;res;0;kind;cut;exp
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

