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
if.1=$y do.
 'dist meth norm data'=: >y
else.
 'dist meth norm data'=: y
end.
 distmx=: dist distancematrix ". norm {:: 'data';'zscoreData data'
 NB. distmx=: dist distancematrix data
 'nr nc' =: $ data
 cmd=. LIBCLUST,' treecluster * c c i i x *x'
 HC=: 0 pick cmd cd dist;meth;nr;nc;distmx;(unDoub data)
 ncl=: 0
 clustid=: nr $ 0
 freedc =: 1 NB. have to stick this in here for .... pointer reasons
 1
)

destroy=: 3 : 0
 if. freedc=1 do. 
 nr freedistmx distmx
 freenodes HC
end.
 freedc=: 0
)

NB. gives cluster labels
cutree=: 3 : 0
 cmd=. LIBCLUST,' cuttree n i x i *i'
 4 pick cmd cd nr;HC;y;clustid
)

uncent=: 4 : 0
 cmd =. LIBCLUST,' farthest_distance d i x'
 0 pick cmd cd y;x
)

cent=: 4 : 0
 cmd =. LIBCLUST,' summed_distances d i x' 
 0 pick cmd cd y;x
)


NB. dumps the tree and splits
dumptree=: 3 : 0
 cmd=. LIBCLUST,' dumpTree i i x *d *i *i'
 nnr =. nr  - 1
 dst =.  memu 0.1 * i.nnr 
 wx =. memu 0 + nnr $ 1-1
 wy =. memu 0 + nnr $ 1-1
'dst wx wy' =. (3,4,5){   cmd cd nnr;HC;dst;wx;wy
 dst;wx,.wy
)


returnClustDx=: 3 : 0
 ncl=: ". ": nclustLogMax 0 pick dumptree a:   NB. weird you need to ".": this
 clustid=:  cutree ncl
)

NB. feed this the tree distances from dumptree
nclustLogMax =:  1+ # - [: maxdx [: diff ^. 

nclustLogMax_z_ =: nclustLogMax_jcluster_

mean=: +/%#
variance=: mean@:*: - *:@mean
zscoreData=:  (-"_1 _ mean) %"_1 _ %:@:variance 
diff=: (] , {:) - {. , ]
maxdx=: [: I. >./ = ]

3 : 0''
if. UNAME-:'Linux' do.
  LIBCLUST=:  '/home/scott/src/jstuff/j-cluster/libhcluster.so'
elseif. UNAME-:'Darwin' do.
  'platform not supported' 13!:8[10
elseif. UNAME-: 'Win' do.
  'platform not supported' 13!:8[10
elseif. do.
  'platform not supported' 13!:8[10
end.
)


unDoub =: 3 : 0
 (15!:14 <'y') +(8*{:$y)*i.{.$y
)

NB. ########################################################## 
NB. two examples of how to use the dll
treetst =: 3 : 0
 ('e';'s') treetst y
:
 'dist meth' =. x
 'nr nc' =. $ y
 cmd=. LIBCLUST,' treecluster * c c i i *x *x'
 0 pick cmd cd dist;meth;nr;nc;(<0);(unDoub y)
)

treetst2 =: 3 : 0
 ('e';'s') treetst2 y
:
 'dist meth' =. x
 'nr nc' =. $ y
 mydist =. dist distancematrix y
 cmd=. LIBCLUST,' treecluster * c c i i x *x'
 out=. 0 pick cmd cd dist;meth;nr;nc;mydist;(unDoub y)
 nr freedistmx mydist
 out
)

NB. dumps the splits and distances for  use with treetst
dumptreeIn=: 4 : 0
cmd=. LIBCLUST,' dumpTree i i x *d *i *i'
 xx =. x - 1 
 dst =.  memu 0.1 * i.xx 
 wx =. memu xx $ 1-1
 wy =. memu xx $ 1-1
'dst wx wy' =. (3,4,5){   cmd cd xx;y;dst;wx;wy
 dst;wx,.wy
)

NB. for treetst
cutreeIn =: 3 : 0
 'nelements tree nclust'=.y
 clustid =. nelements $ 0
 cmd=. LIBCLUST,' cuttree n i x i *i'
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
 'nr nc' =. $ y
 cmd=. LIBCLUST,' distancematrix * i i c *x'
 0 pick cmd cd nr;nc;x;(unDoub y)
)

NB. this frees the allocated distance matrix
freedistmx =: 4 : 0
 cmd=. LIBCLUST,' freedistmx i x x'
 0 pick cmd cd x;y
)

NB. x should be nr
showdists=: 4 : 0
 cmd =. LIBCLUST,' show_dists n i x'
 cmd cd x;y
)

NB. clusterdistance =: 4 : 0
NB.  'nr nc' =. $ y
NB.  wts =. nc $ 2.7 - 1.7
NB.  mask=. masker y
NB.  cmd =. LIBCLUST,' clusterdistance d i i *x *x *d i i *i *i c c i'
NB. )
NB. NB. probably also not needed
NB. somcluster =: 3 : 0
NB.  cmd=. LIBCLUST,' somcluster n x x *d *x *d x x x d x c *d x'
NB. )
NB. NB. potentially useful for comparison to lapack to see what your overhead is like
NB. pcaC =: 3 : 0
NB.  'nr nc' =. $ y
NB.  nn =. (nr <. nc) 
NB.  u =. y
NB.  v =. (nn,nn) $ 1.1 - 1.1
NB.  w =. (nn,nn) $ 1.1 - 1.1
NB.  wts =. ({: $ y) $ 2.7 - 1.7
NB.  cmd =. LIBCLUST,' pca i i i *x *x *x'
NB.  cmd cd nr;nc;(unDoub u);(unDoub v);(unDoub w)
NB.  u;v;w
NB. ) 


NB. NB. henry rich's implementation of  Levenshtein  distance from j-list
NB. levdist=: 4 : 0"1
NB. 'a b'=. (/: #&>)x;y
NB. z=. >: iz =. i.#b
NB. for_j. a do.
NB.    z=. <./\&.(-&iz) (>: <. (j ~: b) + |.!.j_index) z
NB. end.
NB. {:z
NB. ) 
NB. NB. Henry Rich does Jaccard
NB. tanimoto =: (+&#   %/@:-   2 1&*@:(+/@:e.))&~."1 

NB. NB. JP Jacobs
NB. euclidean=: +/&:*:@(-"1)/ 


NB. NB. Further reduction in execution time can be gotten by using integer arithmetic:
NB. NB. d=: +/&:*:@(-"1)/&([: <. (2^20) * ]) 

NB. upTri=: , #~ [: , [: -.@>:/~@i. #

