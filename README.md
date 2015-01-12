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

