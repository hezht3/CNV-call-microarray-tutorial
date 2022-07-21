# ----------------------------- QuantiSNP HowTo ----------------------------- #

cd /users/zhe/GEARS_CNV/quantisnp/

# copy PennCNV tutorial data file
mkdir INPUT
cp /users/zhe/GEARS_CNV/PennCNV-1.0.5/INPUT/sample*.txt INPUT

for f in sample1.txt sample2.txt sample3.txt
do
    awk 'BEGIN{ print "Name" "Chr" "Position" "Log R Ratio" "B Allele Freq"} \
    NR>1 { print $1, $2, $3, $5, $6}' $f \
    > temp.txt; mv temp.txt "$f"
done

mkdir OUTPUT