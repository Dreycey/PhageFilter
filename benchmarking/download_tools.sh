# download and install FVE
git clone https://code.vt.edu/saima5/FastViromeExplorer.git;
cd FastViromeExplorer/;
javac -d bin src/*.java;
cd tools-mac/;
sudo cp kallisto /usr/local/bin/;
cd ../../;

# download and install samtools
git clone https://github.com/samtools/htslib.git;
cd htslib;
git submodule update --init --recursive;
make;
cd ../;
git clone https://github.com/samtools/samtools.git;
cd samtools/;
make;
make install;
sudo cp samtools /usr/local/bin/;
cd ../;

# download and install kraken2 (using conda, since it makes installation easy)
# https://conda.io/projects/conda/en/latest/user-guide/install/macos.html
conda install -c bioconda kraken2;


