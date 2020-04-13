load'tables/csv'
load'../j-cluster.ijs'
raw =. readcsv '../test/usarrests.csv'
data =. > ". each }."1 }.raw
genes=. > ". each readcsv '../test/gene.csv'
rownames =. >0{"1 }.raw
tdrnd=: (([: > [: {. [) + [: ?. 0 $~ ]) ,. ([: > [: {: [) + [: ?. 0 $~ ]
c4=:clust4=: ((50;10)tdrnd 50),((50;_10)tdrnd 50),((100;_90)tdrnd 50),(_90;0)tdrnd 50
NB.clust4 writecsv 'test/4clust.csv'

hr =: ('e';'s') treetst_jcluster_ c4
]dst=.0 pick 200 dumptreeIn_jcluster_ hr

]dst=.0 pick  dumptree_jcluster_ hr
nclustLogMax dst
cutreeIn_jcluster_ 200;hr;4
hr =: ('e';'s') treetst2_jcluster_ c4
dst=.0 pick 200 dumptreeIn_jcluster_ hr
nclustLogMax dst

hr =: ('e';'s') treetst_jcluster_ genes
0 pick 13 dumptreeIn_jcluster_ hr
cutreeIn_jcluster_ 13;hr;4



hc=: conew 'jcluster'
create__hc 'e';'s';0;c4
'dst cmx'=.dumptree__hc ''
nclustLogMax dst
destroy__hc ''




hc=: conew 'jcluster'
create__hc 'e';'s';0;genes
'dst cmx'=.dumptree__hc ''
nclustLogMax dst
destroy__hc ''






hc=: ('e';'s';0;clust4) conew  'jcluster'
cutree__hc 4
'dst cmx'=. dumptree__hc ''
nclustLogMax dst
destroy__hc ''

hc=: ('e';'s';0;clust4) conew 'jcluster'
'dst cmx'=. dumptree__hc ''
cutree__hc 4
destroy__hc''

tsr=: 6!:2 , 7!:2@]    NB. time/space for execution

ndata =. 5e3 800 $?. 4e6$0
hc=: conew 'jcluster'
a=.'e';'s';0;ndata
tsr 'create__hc a'  NB. 8.04126 5.25773e6 for 20k pts, 0.476563 1.32557e6 for 5k
I.1=(cutree__hc 3)
I.0=(cutree__hc 3)
destroy__hc ''

ndata =. 1e4 10 $?. 1e5$0
ndata =. 1e4 10 $ ?. 1e3 $ 0
ndata =: 1e2 10 $ ?. 1e3 $ 0
ndata =: 1e2 10 $ ?. 1e2 $ 0

testmem=: 4 : 0
 for_j. (i.x) do.
  hc=: conew 'jcluster'
  create__hc 'e';'s';0;y
  I.0=cutree__hc 20
  smoutput I.0=cutree__hc 2
  smoutput I.1=cutree__hc 4
  'dst junk'=.dumptree__hc ''
  smoutput ' ' NB. I.0=dst
  destroy__hc ''
 end.
)

testmemNoObj =: 4 : 0
 'nr nc' =. $y
 for_j. (i.x) do.
  tmp2 =. cutreeIn_jcluster_ nr;(('e';'s')treetst_jcluster_ y);2
  tmp4 =. cutreeIn_jcluster_ nr;(('e';'s')treetst_jcluster_ y);4
  smoutput I.0=tmp2
  smoutput I.1=tmp4
  NB.   'dst junk'=.dumptree__hc ''
  smoutput ' ' NB. I.0=dst
 end.
)

testmemNoObj2 =: 4 : 0
 'nr nc' =. $y
 for_j. (i.x) do.
  tmp2 =. cutreeIn_jcluster_ nr;(('e';'s')treetst2_jcluster_ y);2
  tmp4 =. cutreeIn_jcluster_ nr;(('e';'s')treetst2_jcluster_ y);4
  smoutput I.0=tmp2
  smoutput I.1=tmp4
  NB.   'dst junk'=.dumptree__hc ''
  smoutput ' ' NB. I.0=dst
 end.
)

testDistNoObj =: 4 : 0
 'nr nc' =. $ y
 for_j. (i.x) do.
 hr =. ('e';'s') treetst_jcluster_ y
  smoutput I.0 = cutreeIn_jcluster_ (nr;hr;2)
  smoutput I.1 = cutreeIn_jcluster_ (nr;hr;4)
   dist=: 0 pick nr dumptreeIn_jcluster_ hr
NB.  smoutput nclustLogMax 0 pick nr dumptreeIn_jcluster_ hr
 end.
)
1 testDistNoObj genes
hr =: ('e';'s') treetst_jcluster_ genes

distReader=: 'memr x,y,1,8'
hr&distReader 11{ 8 + (16 * i.{.nr)

hr =: ('e';'s') treetst_jcluster_ c4
hr&distReader 3{ 8 + (16 * i.{.nr)
NB. OK, for whatever reason, the shyeah ain't getting in there...
hr =: ('e';'s') treetst_jcluster_ 40}. 55{.c4
hr&distReader 13{ 8 + (16 * i.{.15)

20 testmem ndata
ndata =. 1e4 10 $?. 1e5$0
ndata =. 1e4 10 $ ?. 1e3 $ 0
ndata =: 1e2 10 $ ?. 1e3 $ 0
ndata =: 1e2 10 $ ?. 1e2 $ 0
ndata =: 1e2 10 $ ?. 1e4 $ 0
2 testmem ndata
2 testmemNoObj ndata
2 testmemNoObj2 ndata


dist =. data euclidean data



pdist =. distancematrix data
($ data) showdists pdist

I.0=a=.cutreeIn_jcluster_ 50;(('e';'s')treetst_jcluster_ data);2
I.0=a=.cutreeIn_jcluster_ 50;(('e';'s')treetst2_jcluster_ data);2
I.0=a=.cutreeIn_jcluster_ 50;(('e';'s';'')treetst_jcluster_ data);3
I.0=a=.cutreeIn_jcluster_ 50;(('y';'s')treetst_jcluster_ data);3
I.0=a=.cutreeIn_jcluster_ 50;(('n';'s')treetst_jcluster_ data);3
I.0=a=.cutreeIn_jcluster_ 50;(('o';'s')treetst_jcluster_ data);3

a=. distancematrix_jcluster_ data
50 freedistmx_jcluster_ a
('e';1.1;1.1) calcweights_jcluster_ data 
('e';0.3;0.2) calcweights_jcluster_ data  NB. WTF does this do?





cutreeIn_jcluster_ 50;(('e';'s')treetst_jcluster_ data);10
cutreeIn_jcluster_ 50;(('e';'m')treetst_jcluster_ data);10
cutreeIn_jcluster_ 50;(('e';'a')treetst_jcluster_ data);10
I.0 = a =. cutreeIn_jcluster_ 50;(('e';'c')treetst_jcluster_ data);10 



dat=: > ". each readcsv 'test/gene.csv'
tmp =. ('e';'s')treetst dat
out=: cutree 12;tmp;3
