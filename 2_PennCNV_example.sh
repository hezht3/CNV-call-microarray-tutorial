# ----------------------------- PennCNV user guide ----------------------------- #

cd /users/zhe/GEARS_CNV/PennCNV-1.0.5/example
module load perl
export PATH="/users/zhe/perl5/perlbrew/build/perl-5.14.2-PIC/perl-5.14.2/:$PATH"
export PATH="/users/zhe/GEARS_CNV/PennCNV-1.0.5/:$PATH"

perl runex.pl

# Exercise 1: individual-based calling and write the output to ex1.rawcnv
# Running command <
# detect_cnv.pl \
#    -test \                               # test HMM model to identify CNV
#    -hmm example.hmm \                    # HMM model file
#    -pfb example.pfb \                    # population frequency for B allele file
#    -conf -log ex1.log \                  # calculate confidence for each CNV (experimental feature)
#    -out ex1.rawcnv \                     # specify output root filename
#    father.txt mother.txt offspring.txt   # inputfile
# >
perl runex.pl 1
cat ex1.rawcnv

# Exercise 2: tiro-based calling and write the output to ex2.triocnv
# Running command <
# detect_cnv.pl \
#    -trio \                               # posterior CNV calls for father-mother-offspring trio
#    -hmm example.hmm \                    # HMM model file
#    -pfb example.pfb \                    # population frequency for B allele file
#    -log ex2.log \                        # write notification/warning messages to this file
#    -out ex2.triocnv \                    # specify output root filename
#    -cnv ex1.rawcnv \                     # specify CNV call file for use in family-based CNV calling
#    father.txt mother.txt offspring.txt   # inputfile
# >
perl runex.pl 2
cat ex2.triocnv
cat ex2.log

# Skip exercise 3, no hh550.hg18.gcmodel file

# Exercise 4: validation-based calling for the region between rs8114969 and rs682562
# Running command <
# detect_cnv.pl \
#    -validate \                           # validate copy number at a pre-specified region
#    -hmm example.hmm \                    # HMM model file
#    -pfb example.pfb \                    # population frequency for B allele file
#    -log ex4.log \                        # write notification/warning messages to this file
#    -out ex4.rawcnv \                     # specify output root filename
#    -startsnp rs8114269 \                 # start SNP of a pre-specified region for --validate operation
#    -endsnp rs682562 \                    # end SNP of a pre-specified region for --validate operation
#    -delfreq 0.005 \                      # prior duplication frequency of a pre-specified region for --valiate operation
#    -list inputlist                       # a list file containing path to files to be processed
# >
perl runex.pl 4
cat ex4.rawcnv
cat ex4.log

# Exercise 5: validation-based calling for all CNV regions annotated in 'ccnvr' file
# Running command <
# detect_cnv.pl \
#    -validate \                           # validate copy number at a pre-specified region
#    -hmm example.hmm \                    # HMM model file
#    -pfb example.pfb \                    # population frequency for B allele file
#    -log ex5.log \                        # write notification/warning messages to this file
#    -out ex5.rawcnv \                     # specify output root filename
#    -candlist ccnvr \                     # a file containing all candidate CNV regions to be validated
#    -list inputlist                       # a list file containing path to files to be processed
# >
perl runex.pl 5
cat ex5.rawcnv
cat ex5.log

# Exercise 6: joint CNV calling for a trio and write to ex6.jointcnv
# Running command <
# detect_cnv.pl \
#    -joint \                              # joint CNV calls for trio
#    -hmm example.hmm \                    # HMM model file
#    -pfb example.pfb \                    # population frequency for B allele file
#    -log ex6.log \                        # write notification/warning messages to this file
#    -out ex6.jointcnv \                   # specify output root filename
#    father.txt mother.txt offspring.txt   # inputfile
# >
perl runex.pl 6
cat ex6.jointcnv
cat ex6.log

# Exercise 7: convert CNV call in ex1.rawcnv to BED format and write to ex1.bed
# Running command <
# visualize_cnv.pl \
#    -intype cnv \                         # input file type: cnv (default) or assoc
#    -format bed \                         # out file format (bed, html, beadstudio, wig)
#    -track 'example track' \              # name for UCSC genome browser track (default=inputfile name)
#    -out ex1.bed \                        # output root file name (default=inputfile name)
#    ex1.rawcnv                            # inputfile
# >
perl runex.pl 7
cat ex1.bed

# Exercise 8: convert CNV call in ex3.rawcnv to tab-delimited format and write to ex1.tabcnv
# Running command <
# convert_cnv.pl \
#    -intype penncnv \                     # input file type: cnv (default) or assoc
#    -outtype tab \                        # output file format (default: tab)
#    -output ex1.tabcnv \                  # specify output file name (default: STDOUT)
#    ex1.rawcnv                            # inputfile
# >