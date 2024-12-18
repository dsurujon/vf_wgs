#!/bin/bash

# bcftools
cd ~/tools/
# download the tar: https://www.htslib.org/download/
bzip2 bcftools-1.21.tar.bz2 -d
tar -xf bcftools-1.21.tar 
cd bcftools-1.21
./configure --prefix=/root/bin/bcftools
make
make install

# htslib
cd ~/tools/
# download the tar: https://www.htslib.org/download/
bzip2 htslib-1.21.tar.bz2 -d
tar -xf htslib-1.21.tar 
cd htslib-1.21
./configure --prefix=/root/bin/htslib
make
make install


