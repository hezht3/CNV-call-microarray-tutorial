# ----------------------------- PennCNV user guide ----------------------------- #

cd /users/zhe/GEARS_CNV/PennCNV-1.0.5/example
module load perl
module load conda_R
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
perl runex.pl 8
cat ex1.tabcnv

# Exercise 9: convert CNV call in ex1.newcnv from tab-delimited format to PennCNV format
# Running command <
# convert_cnv.pl \
#    -intype tab \                         # input file type: cnv (default) or assoc
#    -outtype penncnv \                    # specify output file name (default: STDOUT)
#    ex1.tabcnv                            # inputfile
# >
perl runex.pl 9
cat ex1.tabcnv

# Excisise 10: filter CNV calls in ex1.rawcnv and print out only deletions with >=10 SNPs and >=50kb
# Running command <
# filter_cnv.pl \
#    -numsnp 10 \                          # minimum number of SNPs in CNV calls
#    -length 50k \                         # minimum length of CNV calls (suffix of k or m is okay)
#    -type del \                           # type of CNVs (dup or del)
#    ex1.rawcnv                            # inputfile
# >
perl runex.pl 10

# Exercise 11: compare CNV calls between pairs of samples
# Running command <
#    compare_cnv.pl \
#    compdup ex1.rawcnv \                  # CNV calls on duplicated samples in the same file ('dup' operation)
#    -list list.compdup                    # a tab-delimited 2-column file with two identifiers
# >
perl runex.pl 11

# Exercise 12: compare CNV calls between pairs of samples
# Running command <
# compare_cnv.pl \
#    compcall ex1.rawcnv \                 # CNV calls on same sample called by different algorithms in differnet files ('compcall' operation)
#    -list list.compcall \                 # a tab-delimited 2-column file with two identifiers
#    -cnv2 ex2.triocnv                     # callfile 2
# >
perl runex.pl 12

# Exercise 13: generate CNV-based genotype calls)
# Running command <
# infer_snp_allele.pl \
#    -pfb example.pfb \                    # population frequency of B allele file
#    -hmm example.hmm \                    # the HMM file used in PennCNV calling
#    -allcn 221 \                          # copy number of nuclear family (for inherited CNV)
#    -start rs11716390 \                   # start SNP of the CNV
#    -end rs17039742 \                     # end SNP of the CNV
#    -out ex13.geno \                      # the output file (default: STDOUT)
# father.txt mother.txt offspring.txt      # inputfile
# >
perl runex.pl 13
cat ex13.geno

# Exercise 14: validate de novo CNVs and assign P-values
# Running command <
# infer_snp_allele.pl \
#    -pfb example.pfb \                    # population frequency of B allele file
#    -hmm example.hmm \                    # the HMM file used in PennCNV calling
#    -denovocn 1 \                         # copy number of offspring (for de novo CNV)
#    -start rs11716390 \                   # start SNP of the CNV
#    -end rs17039742 \                     # end SNP of the CNV
#    -out ex14.geno \                      # the output file (default: STDOUT)
#    father.txt mother.txt offspring.txt
# >
perl runex.pl 14
cat ex14.geno

# Exercise 15: convert Canary CNV calls to PennCNV format
# Running command <
# convert_cnv.pl \
#    -intype canary \                                # input file format (default: penncnv)
#    -outtype penncnv \                              # output file format (default: tab)
#    -canarydef GenomeWideSNP_6.hg18.cnv_defs \      # file containing CNP locations (for Canary calls)
#    -output ex15.rawcnv \                           # specify output file name (default: STDOUT)
#    example.canary_calls                            # inputfile
# >
perl runex.pl 15
cat ex15.rawcnv

# Exercise 16: generate plots of LRR and BAF values for CNV calls in JPG format (R must be installed in the system for plotting)
# Running command <
# visualize_cnv.pl \
#    -format plot \                        # out file format (bed, html, beadstudio, wig)
#    -signal offspring.txt ex1.rawcnv
# >
perl runex.pl 16