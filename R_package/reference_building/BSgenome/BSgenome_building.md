# BSgenome reference building (not masked !)
## preview
- because of the chromosome names or the version of references of different species, we need to build the necessary BSgenome objects
- the goal of the tutorial is build up BSgenome R packages and install
  
## preparation
- RAM: >=4G if the genome is big, >=2G if small
- software
  - faToTwoBit: convert fasta file to 2bit format
- R package
  - Biostrings
  - BSgenome
- file
  - unmasked and uncompressed `.fa` file (example: [ensembl GRCm38 reference fa file](http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz))

## process
### install R packages
```R
BiocManager::install("Biostrings")
BiocManager::install("BSgenome")
```

### split and convert the reference fa file
- the raw fa file should not in the same dir with the split fa files
- all fa and 2bit files need to be build in the R packages should be in the same dir
```perl
#!/usr/bin/perl                                                                                                                    
# this script is for spliting the reference fa file to each chromosome fa files
# author: laojp
# time: 2023.09.06
# position: SYSUCC bioinformatic platform
# usage
#  enter the result dir and run `perl ....pl [reference file] [output]`
#  then the result will output to the ./

$f = $ARGV[0]; # get the file name
$out = $ARGV[1]; # get the output dir

open (INFILE, "<$f") or die "Can't open: $f $!";

while (<INFILE>) {
        $line = $_; 
        chomp $line;
        if ($line =~ /\>/) { #if > as head
                close OUTFILE;
                my @F = split(/ /,$line);
                my @G = split(/\>/,$F[0]);
                print("$G[1] started \n")
                $line = $F[0];
                $new_file = $out;
                $new_file .= "/";
                $new_file .= $G[1];
                $new_file .= ".fa";
                open (OUTFILE, ">$new_file") or die "Can't open: $new_file $!";
        }
        print OUTFILE "$line\n";
}
close OUTFILE;
```
```bash
# this script is for making twobit files for each fa files 
# author: laojp
# time: 2023.09.06
# position: SYSUCC bioinformatic platform
# usage:
#  fatotwobit_v1.0.sh [input_dir] [output_dir]

# install the software
if [ $(command -v faToTwoBit) ]; then
	echo "faToTwoBit is installed"
else
	mamba install -c bioconda -c conda-forge ucsc-fatotwobit --yes
fi

# scan fa in inputdir and convert to twobit
find ${1} -name "*.fa" | while read id ; do
	echo "${id} started"
	faToTwoBit $(realpath ${id}) ${2}/$(basename ${id} .fa).2bit &
done
```

### create the seed file 
- mainly for reading in the `forgeBSgenomeDataPkg` function
- format: debian control file (DCF)
- content
  - standard description fields
    - Package, Title, Description, Version, Author, Maintainer, License, Suggests
  - not-standard description fields
    - organism, common_name, genome, provider, release_date, source_url, organism_biocview
  - other description fields
    - BSgenomeObjname, seqnames, sirc_seqs, mseqnames, seqs_srcdir, seqfile_name, seqfiles_prefix, seqfiles_suffix, ondisk_seq_format
```txt
# the content is the example for the seed file, name: BSgenome.Mmusculus.ENSEMBL.V102-seed

Package: BSgenome.Mmusculus.ENSEMBL.V102
Title: BSgenome for ensembl GRCm38.p6 (V102)
Description: the BSgenome reference packages for ensembl GRCm38.p6 (version 102)
Version: 1.0
Author: champeil
Maintainer: champeil
License: champeil in SYSUCC bioinformatic platform
organism: Mus musculus
common_name: Mouse
genome: GRCm38.p6
provider: ENSEMBL
release_date: 2020-11
source_url: http://ftp.ensembl.org/pub/release-102/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.primary_assembly.fa.gz
organism_biocview: Mus_musculus
BSgenomeObjname: Mmusculus
seqnames: gsub(".fa","",list.files("/home/laojp/database/mouse/ensembl/fasta/GRCm38/ensembl_102/bsgenome/",pattern=".fa"))
seqs_srcdir: /home/laojp/database/mouse/ensembl/fasta/GRCm38/ensembl_102/bsgenome/
ondisk_seq_format: 2bit
```

### read seed and build source tree
```R
# this script is for building source tree
# author: laojp
# time: 2023.09.06
# position: SYSUCC bioinformatic platform

library(Biostrings)
library(BSgenome)

forgeBSgenomeDataPkg("/home/laojp/database/mouse/ensembl/fasta/GRCm38/ensembl_102/bsgenome/BSgenome.Mmusculus.ENSEMBL.V102-seed")
```
### read source tree and build-check-install packages
```bash
R CMD build <pkg_dir: # path to the source tree offf the package>
R CMD check <tarball: # path to the tarball builded by R CMD>
R CMD INSTALL <tarball: # path to the tarball builded by R CMD>>
```
## summary and the code
- download the `make_BSgenome.sh`, `build_source_tree.R`, `split_reference.pl` scripts in the same dir
- prepare the necessary input for `make_BSgenome.sh` and run
- the chinese version tutorial is in my obsidian md file, contact if you want

## a big pie in the sky
- will update the masked BSgenome
- will update the TxDb version 

## citation
if you like my tutorial, please cite:
- Zhang M, Wen H, Liang M, Qin Y, Zeng Q, Luo D, Zhong X, Li S. Diagnostic Value of Sylvian Fissure Hyperechogenicity in Fetal SAH. AJNR Am J Neuroradiol. 2022 Apr;43(4):627-632. doi: 10.3174/ajnr.A7449. Epub 2022 Mar 10. PMID: 35272984; PMCID: PMC8993207.

## reference
- https://github.com/mevers/build_custom_BSgenome_TxDb
- https://www.bioconductor.org/packages/release/bioc/vignettes/BSgenome/inst/doc/BSgenomeForge.pdf
- https://blog.csdn.net/hs6605015/article/details/119819558
- https://www.jianshu.com/p/6073a75870d5










