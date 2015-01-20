load'tables/csv'
load'j-cluster.ijs'
raw =. readcsv 'test/usarrests.csv'
data =. > ". each }."1 }.raw
rownames =. >0{"1 }.raw


hc=: conew 'jcluster'
create__hc 'e';'s';data
cutree__hc 20
destroy__hc ''


ndata =. 2e4 10 $?. 2e5$0
hc=: conew 'jcluster'
create__hc 'e';'s';ndata
I.1=(cutree__hc 3)
I.0=(cutree__hc 3)
destroy__hc ''

ndata =. 1e4 10 $?. 1e5$0
testmem=: 4 : 0
for_j. (i.x) do.
hc=: conew 'jcluster'
create__hc 'e';'s';y
I.0=cutree__hc 20
I.0=cutree__hc 2
smoutput I.1=cutree__hc 100
destroy__hc ''
end.
)

20 testmem ndata
ndata =. 2e4 10 $?. 2e5$0
20 testmem ndata

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
