#!/bin/bash
make

# build inhouse components
g++ -DHAVE_CONFIG_H -I.  -O3   -g -O2 -MT art_illumina_src/cvgdist.o -MD -MP -c -o art_illumina_src/cvgdist.o art_illumina_src/cvgdist.cpp
# build executable using inhouse components
g++  -g -O2   -o art_illumina art_illumina_src/art_illumina.o art_illumina_src/art_qual_scale.o art_illumina_src/empdist.o art_illumina_src/readSeqFile.o art_illumina_src/seqRead.o art_illumina_src/samRead.o art_illumina_src/cvgdist.o  -lgsl -lgslcblas -lm

make

# test me!
#./art_illumina -sam -p -l 150 -ss HS25 -m 500 -s 10 -i genome.fa -f 10 -o reads -cp test.cfa
