### AR Bank scripts
#### Bin, lib, and script dependencies:
- Biopython
- Bandage
- BBTools
- BLAST+
- c-SSTAR
- Kraken
- mlst
- Prokka
- QUAST
- QualAssessCleanSeqs.bash
- QualAssessRawSeqs.bash
- RNAmmer
- SPAdes
- summarize_kraken-report.sh
- Trimmomatic

#### Ref sequences and database dependencies:
- PhiX genome (FastA format)
- Illumina adapter sequences (FastA format)
- ARG-ANNOT AR database (SRST2-formatted FastA)
- plasmids database (FastA format)
- a single full-length 16S rRNA gene sequence (FastA format)

#### On SGE cluster:

`module load Bandage/0.7.1 bbmap trimmomatic/0.35 kraken/0.10.5 SPAdes/3.6.2 quast/2.3 ncbi-blast+/2.2.30 prokka/1.8 rnammer/1.2 Python/2.7.3`

#### Example local installation:

    brew tap homebrew/science && brew update
    brew install bbmap blast bowtie2 kraken mlst prokka quast rnammer spades trimmomatic
    brew tap tseemann/homebrew-bioinformatics-linux && brew update
    brew install bandage
    pip install -U pip && pip install biopython
    cd $HOME
    git clone https://github.com/chrisgulvik/AR_Bank_scripts.git
    echo 'export PATH="$PATH:$HOME/AR_Bank_scripts"' >> $HOME/.bash_profile
    gunzip ~/AR_Bank_scripts/DBs/*.gz
    echo 'export AR_BANK_DBs="$HOME/AR_Bank_scripts/DBs"' >> $HOME/.bash_profile
    git clone https://github.com/chrisgulvik/c-SSTAR.git
    echo 'export PATH="$PATH:$HOME/c-SSTAR"' >> $HOME/.bash_profile
    git clone https://github.com/chrisgulvik/genomics_scripts.git
    echo 'export PATH="$PATH:$HOME/genomics_scripts"' >> $HOME/.bash_profile
    git clone https://github.com/chrisgulvik/summarize_kraken_data.git
    echo 'export PATH="$PATH:$HOME/summarize_kraken_data"' >> $HOME/.bash_profile
    source ~/.bash_profile
