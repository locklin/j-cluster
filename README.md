j-cluster
=========

Heirarchical clustering for J, taken from this old biosciences library:

http://bonsai.hgc.jp/~mdehoon/software/cluster/software.htm#ctv

Methods beyond single linkage are presently untested.

The pointer to pointer stuff is readable, but kind of silly.

For now this is a nice template to develop new Hierarchical clustering
distance metrics.

The mask piece had to be removed due to a peculiarity with J's FFI which I was 
unable to overcome. It was a bloody mess and makes no mathematical sense 
anyway. 

Having the flexibility to generate curated distance metrics from J might
be helpful at some point, so I'll leave the logic for that for now.


