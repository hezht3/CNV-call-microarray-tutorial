# ----------------------------- Concordant call anneal ----------------------------- #

cd /users/zhe/GEARS_CNV/PennCNV-1.0.5/
module load perl
module load conda_R
export PATH="/users/zhe/perl5/perlbrew/build/perl-5.14.2-PIC/perl-5.14.2/:$PATH"
export PATH="/users/zhe/GEARS_CNV/PennCNV-1.0.5/:$PATH"

# convert concordant calls to PennCNV format
convert_cnv.pl \
    --intype tab \
    --outtype penncnv \
    ./OUTPUT/sampleall.concordant.rawcnv \
    > ./OUTPUT/sampleall.concordant.penncnv
rm ./OUTPUT/sampleall.concordant.rawcnv

# anneal
clean_cnv.pl \
    combineseg \
    --signalfile INPUT/sample_pfb.pfb \
    ./OUTPUT/sampleall.concordant.penncnv \
    > ./OUTPUT/sampleall.concordant.anneal.penncnv