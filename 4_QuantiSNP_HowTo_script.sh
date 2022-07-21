#!/bin/csh

cd /users/zhe/GEARS_CNV/quantisnp
#
# EXAMPLE SHELL SCRIPT FOR RUNNING QUANTISNP 2
#

# define program variables
set EMITERS 		= 10									               		# number of EM iterations to use during training
set LSETTING 		= 2000000								            		# characteristic CNV length parameter
set GCDIR 			= /users/zhe/GEARS_CNV/quantisnp/gc_content/b35/			# set path to GC data files (contents of gc_data.zip)
set SUBSAMPLELEVEL  = 1 							            				# number of sub-samples to use
set PARAMSFILE 		= ./config/params.dat								     	# path to parameters file
set LEVELSFILE 		= ./config/levels.dat								    	# path to levels file
set CHRRANGE 		= 1:23											        	# path to parameters file
set CHRX 			= 23											        	# which chromosome is X?

set OUTDIR 			= /users/zhe/GEARS_CNV/quantisnp/OUTPUT/		         	# output directory

set SAMPLEID 		= 99HI0697A										        	# sample name
# set GENDER 			= female									            # sample gender
set INFILE 			= /users/zhe/GEARS_CNV/quantisnp/INPUT/sample1.txt	        # input data file


# set path to MCR Run-Time Libraries
set MCRROOT		= /users/zhe/GEARS_CNV/MATLAB/MATLAB_Compiler_Runtime/v79

mkdir -p $OUTDIR

./linux64/run_quantisnp2.sh \
    $MCRROOT \
    --chr $CHRRANGE \
    --outdir $OUTDIR \
    --sampleid $SAMPLEID \
    #--gender $GENDER \
    --emiters $EMITERS \
    --lsetting $LSETTING \
    --gcdir $GCDIR \
    --plot \
    --genotype \
    --config $PARAMSFILE \
    --levels $LEVELSFILE \
    --input-files $INFILE \
    --chrX $CHRX \
    --doXcorrect

echo "EXAMPLE END"
