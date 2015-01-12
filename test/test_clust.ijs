load'tables/csv'
load'j-cluster.ijs'
raw =. readcsv 'test/usarrests.csv'
data =. > ". each }."1 }.raw
rownames =. >0{"1 }.raw
dist =. data euclidean data

I.0=a=.cutree 50;(('e';'s')treetst data);2
I.0=a=.cutree 50;(('e';'s')treetst data);3
I.0=a=.cutree 50;(('y';'s')treetst data);3
I.0=a=.cutree 50;(('n';'s')treetst data);3
I.0=a=.cutree 50;(('o';'s')treetst data);3

distancematrix data

('e';1.1;1.1) calcweights data 
('e';0.3;0.2) calcweights data  NB. WTF does this do?





cutree 50;(('e';'s')treetst data);10
cutree 50;(('e';'m')treetst data);10
cutree 50;(('e';'a')treetst data);10
I.0 = a =. cuttree 50;(('e';'c')treetst data);10 dumps core

cutree 50;(('e';'a')treetst data);3
cutree 50;(('e';'a')treetst data);3


dat=: > ". each readcsv 'test/gene.csv'
tmp =. ('e';'s')treetst dat
out=: cutree 12;tmp;3
