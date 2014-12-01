load'tables/csv'
load'j-cluster.ijs'
dat=: > ". each readcsv 'test/gene.csv'
(1;0) euclid2 dat
(1;0)mymetric dat



raw =. readcsv 'test/usarrests.csv'
data =. > ". each }."1 }.raw
rownames =. >0{"1 }.raw
dist =. data euclidean data

'e' distancematrix (] ; 1 #~ [: }. $) data
tmp =. treetst2 data
out=: cuttree 50;tmp;10
cuttree 50;(('e';'m')treetst2 data);10
cuttree 50;(('e';'a')treetst2 data);10
NB. cuttree 50;(('e';'c')treetst2 data);10 dumps core

out=:cuttree 50;(('e';'s')treetst2 data);3
out=: cuttree 50;(('e';'m')treetst2 data);3
cuttree 50;(('e';'a')treetst2 data);3


dat=: > ". each readcsv 'test/gene.csv'
tmp =. ('e';'s')treetst2 dat
out=: cuttree 12;tmp;3
