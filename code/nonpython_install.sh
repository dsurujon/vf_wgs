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


# spades
cd ~/tools
wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz
tar -xzf SPAdes-4.0.0-Linux.tar.gz
export PATH=$PATH:/root/tools/SPAdes-4.0.0-Linux/bin/

# seqtk
cd ~/tools
git clone https://github.com/lh3/seqtk.git
cd seqtk
make
cp seqtk /root/bin/