load'tables/csv'
load'j-cluster.ijs'
loc=. 3 : '> (4!:4 <''y'') { 4!:3 $0'  
PATH=: getpath_j_ loc''
raw =. readcsv PATH,'./usarrests.csv'
data =. > ". each }."1 }.raw
rownames =. >0{"1 }.raw

Note 'To run j-cluster tests:'
  load 'j-cluster.ijs'
  load 'test/test_clust.ijs'
)

tdrnd=: (([: > [: {. [) + [: ?. 0 $~ ]) ,. ([: > [: {: [) + [: ?. 0 $~ ]

clust4=: ((50;10)tdrnd 50),((50;_10)tdrnd 50),((100;_90)tdrnd 50),(_90;0)tdrnd 50

hc=: conew 'jcluster'
create__hc 'e';'s';0;clust4
'dst ab'=. dumptree__hc ''
nclustLogMax dst

destroy__hc''
testcutoff =: 3 : 0

)

testc=: 3 : 0
hc=: conew 'jcluster'
create__hc 'e';'s';0;data
  assert. 32=I.0=cutree__hc 2
  assert. 0 2 3 4 5 6 7 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 31 33 34 35 36 37 38 39 40 41 42 43 44 45 46 47 48 49 =I.3=cutree__hc 4
 'dst cmx'=. dumptree__hc ''  
  assert. 371.1= {:dst
  assert. 32 _48={:cmx
  destroy__tree ''
)


smoutput testc ''
smoutput testcutoff ''




