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

# tRNAscan-SE
cd ~/tools
wget https://github.com/UCSC-LoweLab/tRNAscan-SE/archive/refs/tags/v2.0.11.tar.gz
mv v2.0.11.tar.gz tRNAscan-SE.tar.gz
tar -xzf tRNAscan-SE.tar.gz
cd tRNAscan-SE-2.0.11 
./configure --prefix /root/bin 
make 
make install

# Infernal 1.1.4
cd ~/tools
wget http://eddylab.org/infernal/infernal-1.1.4-linux-intel-gcc.tar.gz
tar -xzf infernal-1.1.4-linux-intel-gcc.tar.gz
cp infernal-1.1.4-linux-intel-gcc/binaries/* /root/bin/
#Aragorn
curl --location --remote-name http://www.trna.se/ARAGORN/Downloads/aragorn1.2.41.c
gcc -O3 -ffast-math -finline-functions -o aragorn aragorn1.2.41.c
mv aragorn /root/bin/
#Diamond
curl --location --remote-name https://github.com/bbuchfink/diamond/releases/download/v2.1.10/diamond-linux64.tar.gz
tar xzf diamond-linux64.tar.gz
rm diamond-linux64.tar.gz
mv diamond /root/bin/

#mafft
wget https://mafft.cbrc.jp/alignment/software/mafft-7.525-without-extensions-src.tgz
gunzip -cd mafft-7.525-without-extensions-src.tgz | tar xfv -
cd mafft-7.525-without-extensions/core
make clean
make
su
make install
cp mafft /root/bin/
cd ../..
rm mafft-7.525-without-extensions-src.tgz