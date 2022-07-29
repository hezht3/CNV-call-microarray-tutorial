# ----------------------------- QuantiSNP HowTo ----------------------------- #

cd /users/zhe/GEARS_CNV/quantisnp/

module load perl
module load conda_R
export PATH="/users/zhe/perl5/perlbrew/build/perl-5.14.2-PIC/perl-5.14.2/:$PATH"
export PATH="/users/zhe/GEARS_CNV/PennCNV-1.0.5/:$PATH"

# copy PennCNV tutorial data file
mkdir INPUT
cp /users/zhe/GEARS_CNV/PennCNV-1.0.5/INPUT/sample*.txt INPUT

for f in sample1.txt sample2.txt sample3.txt
do
    awk 'BEGIN{ print "Name" "Chr" "Position" "Log R Ratio" "B Allele Freq"} \
    NR>1 { print $1, $2, $3, $5, $6}' $f \
    > temp.txt; mv temp.txt "$f"
done

# CNV call
mkdir OUTPUT

chmod u+x ./script/4_QuantiSNP_HowTo_script.sh
./script/4_QuantiSNP_HowTo_script.sh

# convert output to PennCNV format
# tab file col: chr, start pos, end pos, cn, sample ID, start SNP, end SNP, conf, numsnp
mkdir ./OUTPUT/convert
awk 'NR>1 {print $2 "\t" $3 "\t" $4 "\t" $9 "\t" $1 "\t" $5 "\t" $6 "\t" $10 "\t" $8}' ./OUTPUT/99HI0697A.cnv > ./OUTPUT/convert/99HI0697A.cnv.penncnv
convert_cnv.pl \
    --intype tab \
    --outtype penncnv \
    ./OUTPUT/convert/99HI0697A.cnv.penncnv \
    > ./OUTPUT/convert/99HI0697A.penncnv

awk 'NR>1 {print $2 "\t" $3 "\t" $4 "\t" $9 "\t" $1 "\t" $5 "\t" $6 "\t" $10 "\t" $8}' ./OUTPUT/99HI0698C.cnv > ./OUTPUT/convert/99HI0698C.cnv.penncnv
convert_cnv.pl \
    --intype tab \
    --outtype penncnv \
    ./OUTPUT/convert/99HI0698C.cnv.penncnv \
    > ./OUTPUT/convert/99HI0698C.penncnv

awk 'NR>1 {print $2 "\t" $3 "\t" $4 "\t" $9 "\t" $1 "\t" $5 "\t" $6 "\t" $10 "\t" $8}' ./OUTPUT/99HI0700A.cnv > ./OUTPUT/convert/99HI0700A.cnv.penncnv
convert_cnv.pl \
    --intype tab \
    --outtype penncnv \
    ./OUTPUT/convert/99HI0700A.cnv.penncnv \
    > ./OUTPUT/convert/99HI0700A.penncnv

# anneal
clean_cnv.pl \
    combineseg \
    --signalfile /users/zhe/GEARS_CNV/PennCNV-1.0.5/INPUT/sample1.txt \
    ./OUTPUT/convert/99HI0697A.penncnv \
    > ./OUTPUT/convert/99HI0697A.adjust.penncnv

clean_cnv.pl \
    combineseg \
    --signalfile /users/zhe/GEARS_CNV/PennCNV-1.0.5/INPUT/sample2.txt \
    ./OUTPUT/convert/99HI0698C.penncnv \
    > ./OUTPUT/convert/99HI0698C.adjust.penncnv

clean_cnv.pl \
    combineseg \
    --signalfile /users/zhe/GEARS_CNV/PennCNV-1.0.5/INPUT/sample2.txt \
    ./OUTPUT/convert/99HI0700A.penncnv \
    > ./OUTPUT/convert/99HI0700A.adjust.penncnv