# ----------------------------- QuantiSNP HowTo ----------------------------- #

cd /users/zhe/GEARS_CNV/quantisnp/

# export MCRROOT="/users/zhe/GEARS_CNV/MATLAB/MATLAB_Compiler_Runtime/v79"

# export LD_LIBRARY_PATH="${LD_LIBRARY_PATH:+${LD_LIBRARY_PATH}:}\
#     /users/zhe/GEARS_CNV/MATLAB/MATLAB_Compiler_Runtime/v79/runtime/glnxa64:\
#     /users/zhe/GEARS_CNV/MATLAB/MATLAB_Compiler_Runtime/v79/bin/glnxa64:\
#     /users/zhe/GEARS_CNV/MATLAB/MATLAB_Compiler_Runtime/v79/sys/os/glnxa64:\
#     /users/zhe/GEARS_CNV/MATLAB/MATLAB_Compiler_Runtime/v79/extern/lib/glnxa64:\
#     /users/zhe/GEARS_CNV/MATLAB/MATLAB_Compiler_Runtime/v79/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:\
#     /users/zhe/GEARS_CNV/MATLAB/MATLAB_Compiler_Runtime/v79/sys/java/jre/glnxa64/jre/lib/amd64/server:\
#     /users/zhe/GEARS_CNV/MATLAB/MATLAB_Compiler_Runtime/v79/sys/java/jre/glnxa64/jre/lib/amd64"

# export XAPPLRESDIR="/users/zhe/GEARS_CNV/MATLAB/MATLAB_Compiler_Runtime/v79/X11/app-defaults/"

# copy PennCNV tutorial data file
mkdir INPUT
cp /users/zhe/GEARS_CNV/PennCNV-1.0.5/INPUT/sample*.txt INPUT

mkdir OUTPUT

chmod u+x ./linux64/run_quantisnp2_modify.sh

./linux64/run_quantisnp2_modify.sh \   # modify `MCRROOT` in `run_quantisnp2.sh`
    --outdir ./OUTPUT/ \
    --levels ./config/levels.dat \
    --config ./config/params.dat \
    --sampleid 99HI0697A \
    --input-files ./INPUT/sample1.txt