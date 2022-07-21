# ----------------------------- Installation ----------------------------- #
qrsh -l mem_free=20G,h_vmem=20G,h_fsize=20G

# ----------------------------- PennCNV ----------------------------- #
mkdir GEARS_CNV
cd GEARS_CNV
wget https://github.com/WGLab/PennCNV/archive/v1.0.5.tar.gz
tar xvfz v1.0.5.tar.gz

module load perl

# compile Perl 5.14.2 using perlbrew
\curl -L https://install.perlbrew.pl | bash   # install perlbrew
perlbrew install perl-5.14.2 --as perl-5.14.2-PIC -Accflags=-fPIC
perlbrew switch perl-5.14.2-PIC

# compile PennCNV
cd PennCNV-1.0.5/kext/
make

# ----------------------------- iPattern ----------------------------- #
cd /users/zhe/GEARS_CNV/

wget https://www.tcag.ca/documents/tools/ipn_0.582.tar.gz
tar xvfz ipn_0.582.tar.gz

export IPNBASE='/users/zhe/GEARS_CNV/ipn_0.582'
PYTHONPATH=$PYTHONPATH:'$IPNBASE/ipnlib'
export PATH="$IPNBASE/preprocess/affy"
export PATH="$IPNBASE/preprocess/ilmn"
export PATH="$IPNBASE/ipn"

cd /users/zhe/GEARS_CNV/ipn_0.582/ipnlib/
ln -sf ipn_pbs_qsub.py ipn_qsub.py

# ----------------------------- QuantiSNP ----------------------------- #
# reference: https://github.com/hernanmd/QuantiSNP2
cd /users/zhe/GEARS_CNV/

# install MCR
wget ftp://ftp.stats.ox.ac.uk/pub/yau/quantisnp2/mcr/MCRinstaller64.run
chmod 755 MCRinstaller64.run
./MCRinstaller64.run  --noexec --keep --nox11 --target Install
cd Install
./MCRInstaller.bin -is:extract
cd istemp43633202141501
chmod 755 jre1.5.0-linux-amd64.bin
./jre1.5.0-linux-amd64.bin
./bin/java -jar setup.jar -console

# install quantisnp
cd ../..
wget ftp://ftp.stats.ox.ac.uk/pub/yau/quantisnp2/executables/09042010/install_quantisnp
chmod 755 install_quantisnp
./install_quantisnp

# download local GC content file
cd quantisnp/
mkdir gc_content
cd gc_content/

wget ftp://ftp.stats.ox.ac.uk/pub/yau/quantisnp2/download/b35.tar.gz
wget ftp://ftp.stats.ox.ac.uk/pub/yau/quantisnp2/download/b36.tar.gz
wget ftp://ftp.stats.ox.ac.uk/pub/yau/quantisnp2/download/b37.tar.gz

ls b*.tar.gz |xargs -n1 tar -xzf