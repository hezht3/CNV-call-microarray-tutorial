# ----------------------------- PennCNV user guide ----------------------------- #

cd /users/zhe/GEARS_CNV/PennCNV-1.0.5/
module load perl
module load conda_R
export PATH="/users/zhe/perl5/perlbrew/build/perl-5.14.2-PIC/perl-5.14.2/:$PATH"
export PATH="/users/zhe/GEARS_CNV/PennCNV-1.0.5/:$PATH"

# download tutorial data file
mkdir INPUT
cd INPUT/
wget --no-check-certificate http://www.openbioinformatics.org/penncnv/download/tutorial_text.zip
unzip tutorial_text.zip

# split data file for each individual
cut -f 1-6 tutorial.txt > sample1.txt
cut -f 1-3,7-9 tutorial.txt > sample2.txt
cut -f 1-3,10-12 tutorial.txt > sample3.txt

kcolumn.pl tutorial.txt split 3 -heading 3 -tab -out sample   # alternative method, write to sample.split*

# CNV call
cd ..
mkdir LOG/
mkdir OUTPUT/

compile_pfb.pl \
    INPUT/sample*.txt \
    -output INPUT/sample_pfb.pfb   # compile population B allele frequency using samples

detect_cnv.pl \
    -test \
    -hmm lib/hh550.hmm \
    -pfb INPUT/sample_pfb.pfb \
    INPUT/sample*.txt \
    -log LOG/sampleall.log \
    -out OUTPUT/sampleall.rawcnv

cat OUTPUT/sampleall.rawcnv
wc -l OUTPUT/sampleall.rawcnv
cat LOG/sampleall.log

### GC model adjustment
gzip -d gc_file/hg18.gc5Base.txt.gz
cal_gc_snp.pl \
    gc_file/hg18.gc5Base.txt INPUT/sample*.txt \
    -output INPUT/sample_gc.gcmodel   # generate GC model file

detect_cnv.pl \
    -test \
    -hmm lib/hh550.hmm \
    -pfb INPUT/sample_pfb.pfb \
    INPUT/sample*.txt \
    -log LOG/sampleall.adjusted.log \
    -out OUTPUT/sampleall.adjusted.rawcnv \
    -gcmodel INPUT/sample_gc.gcmodel

cat OUTPUT/sampleall.adjusted.rawcnv
wc -l OUTPUT/sampleall.adjusted.rawcnv
cat LOG/sampleall.adjusted.log

# anneal
clean_cnv.pl \
    combineseg \
    --signalfile ./INPUT/sample1.txt \
    ./OUTPUT/sampleall.adjusted.rawcnv \
    > ./OUTPUT/sampleall.adjusted.anneal.rawcnv