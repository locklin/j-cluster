j-cluster
=========

Taken from this old biosciences library:

http://bonsai.hgc.jp/~mdehoon/software/cluster/software.htm#ctv

Eventually will change the mask variables to logical (or at least char) 
arrays as they should be.

The pointer to pointer crap is readable, but kind of silly.
For now this is a nice template to develop new Hierarchical clustering
metrics.

I added 3 so far; Chebyshev, Cosine and Angle.

I guess having the flexibility to generate distance metrics from J might
be helpful at some point, so I'll leave the logic for that for now.

One thing which should be noted: the actual distance calculations in the original
library are wrong. Euclidean is anyway. However, it is mathematically self
consistent, so as long as you're not using the absolute distances for something
it should work just fine for clustering.

