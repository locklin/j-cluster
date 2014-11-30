load'tables/csv'
load'j-cluster.ijs'
raw =. readcsv 'test/usarrests.csv'
data =. > ". each }."1 }.raw
rownames =. >0{"1 }.raw
dist =. data euclidean data

'e' distancematrix (] ; 1 #~ [: }. $) data
tmp =. treetst data
out=: cuttree 50;tmp;10
out

