load'tables/csv'
load'j-cluster.ijs'
raw =. readcsv 'test/usarrests.csv'
data =. > ". each }."1 }.raw
rownames =. >0{"1 }.raw
dist =. data euclidean data
I.0=a=.cutree 50;(('e';'s')treetst data);2


b=. 0.1 + (3 3) $ 0 1 2 0 3 8 1 2 3
(b euclidean b)^0.5
cutree 3;(('e';'s')treetst2  b);2
cutree 3;(('e';'s')treetst  b);2
('e';'s')treetst  0.1+b

dat=: > ". each readcsv 'test/gene.csv'
cutree 3;(('e';'s')treetst2  dat);2
(1;0) euclid2 dat
(1;0)mymetric dat



raw =. readcsv 'test/usarrests.csv'
data =. > ". each }."1 }.raw
rownames =. >0{"1 }.raw
dist =. data euclidean data
cutree 50;(('e';'s')treetst2 data);5
cutree 50;(('e';'a')treetst2 data);2
NB. cuttree 50;(('e';'s')treetst data);3

tmp =. treetst2 data
out=: cuttree 50;tmp;10
cuttree 50;(('e';'s')treetst2 data);10
cuttree 50;(('e';'m')treetst2 data);10
cuttree 50;(('e';'a')treetst2 data);10
NB. cuttree 50;(('e';'c')treetst2 data);10 dumps core

out=:cuttree 50;(('e';'s')treetst2 data);3
out=: cuttree 50;(('e';'m')treetst2 data);3
cuttree 50;(('e';'a')treetst2 data);3
cuttree 50;(('e';'a')treetst data);3


dat=: > ". each readcsv 'test/gene.csv'
tmp =. ('e';'s')treetst2 dat
out=: cuttree 12;tmp;3
