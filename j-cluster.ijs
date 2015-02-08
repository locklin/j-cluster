coclass 'jcluster'

NB. 'metric;clustmeth;normalize;data'
NB. available distance metrics:
NB. case 'e': euclid;
NB. case 'b': cityblock;
NB. case 'c': correlation;
NB. case 'a': acorrelation;
NB. case 'u': ucorrelation;
NB. case 'x': uacorrelation;
NB. case 's': spearman;
NB. case 'k': kendall;
NB. case 'o': angle; 
NB. case 'n': cosine; 
NB. case 'y': chebyshev;
NB.
NB. available cluster methods:
NB. 's': pairwise single-linkage clustering
NB. 'm': pairwise maximum- (or complete-) linkage clustering
NB. 'a': pairwise average-linkage clustering
NB. 'c': pairwise centroid-linkage clustering
NB. 
NB. normalize; 0 for none, 1 for variance normalized
NB. 
NB. data; a rank 2 array of doubles
create=: 3 : 0
 'dist meth norm data'=: y
 distmx=: dist distancematrix ". norm {:: 'data';'zscoreData data'
 wt =: ({: $ y) $ 2.7 - 1.7
 mask =.  masker data
 'nr nc' =: $ data
 cmd=. LIBCLUST,' treecluster * x x *x *x *d x c c x'
 HC=: 0 pick cmd cd nr;nc;(unDoub data); (unInt mask) ;wt;0;dist;meth;distmx
 1
)

destroy=: 3 : 0
 nr freedistmx distmx
 freenodes HC
 codestroy ''
)

NB. gives clusters
cutree=: 3 : 0
 clustid =. nr $ 0
 cmd=. LIBCLUST,' cuttree n x x x *i'
 4 pick cmd cd nr;HC;y;clustid
)

uncentral=: 3 : 0
 distmx&uncent"0 i.nr
)

central=: 3 : 0
 distmx&cent"0 i.nr
)

uncent=: 4 : 0
 cmd =. LIBCLUST,' farthest_distance d i x'
 0 pick cmd cd y;x
)

cent=: 3 : 0
 cmd =. LIBCLUST,' summed_distance d i x' 
 0 pick cmd cd y;x
)

NB. dumps the tree and splits
dumptree=: 3 : 0
 cmd=. LIBCLUST,' dumpTree i i x *i *i *d'
 dst =.  (nr-1) $ 1.7 - 1.7
 wx =. (nr-1) $ 1-1
 wy =. (nr-1) $ 1-1
 cmd cd nr;HC;wx;wy;dst
 dst;wx,.wy
)

NB. feed this the tree distances from dumptree
nclustLogMax=: # - [: maxdx [: diff ^.
nclustLogMax_z_ =: nclustLogMax_jcluster_



clustdist=: 3 : 0
 cmd=. LIBCLUST,' clusterdistance d x x x x *d i i *i *i c c x'
)

mean=: +/%#
variance=: mean@:*: - *:@mean
zscoreData=:  (-"_1 _ mean) %"_1 _ %:@:variance 
diff=: (] , {:) - {. , ]
maxdx=: [: I. >./ = ]

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


treetst2 =: 3 : 0
 ('e';'s') treetst y
:
 'dist meth' =. x
 wt =. ({: $ y) $ 2.7 - 1.7
 mask =.  masker y
 'nr nc' =. $ y
 mydist =. distancematrix y
 cmd=. LIBCLUST,' treecluster * x x *x *x *d x c c x'
 out=. 0 pick cmd cd nr;nc;(unDoub y); (unInt mask) ;wt;0;dist;meth;mydist
 nr freedistmx mydist
 out
)



cutreeIn =: 3 : 0
 'nelements tree nclust'=.y
 clustid =. nelements $ 0
 cmd=. LIBCLUST,' cuttree n x x x *i'
 4 pick cmd cd nelements;tree;nclust;clustid
)

NB. frees the malloc in the hclust piece
freenodes=: 3 : 0
 cmd=. LIBCLUST,' freeNodes i x'
 0 pick cmd cd y
)

NB. initializes distance matrix
distancematrix=: 3 : 0
 'e' distancematrix y
:
 wts =. ({: $ y) $ 2.7 - 1.7
 'nr nc' =. $ y
 mask =. masker y
 cmd=. LIBCLUST,' distancematrix * x x *x *x *d c x'
 0 pick cmd cd nr;nc;(unDoub y);(unInt mask);wts;x;0
)

NB. this frees the allocated distance matrix
freedistmx =: 4 : 0
 cmd=. LIBCLUST,' freedistmx i x x'
 0 pick cmd cd x;y
)

showdists=: 4 : 0
'nr nc'=.x
cmd =. LIBCLUST,' show_dists n i x'
cmd cd nr; y
)

NB. probably also not needed
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

