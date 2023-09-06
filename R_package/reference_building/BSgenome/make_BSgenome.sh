# this script is for making BSgenome and install in R
# author: laojp
# time: 2023.09.06
# position: SYSUCC bioinformatic platform
# usage: 
#  make_BSgenome.sh [seed_file] [fa_file] [output_dir]
#   [fa]: all_fa uncompressed .fa file
#  place all file in the same dir
#  the fa_file can not in the same dir with output_dir

script_dir=$(readlink -f "$0" | xargs dirname)

echo "---------- first split of the reference fa ----------"
perl ${script_dir}/split_reference.pl ${2} ${3}

echo "---------- second make the twobit file for all fa file ----------"
if [ $(command -v mamba) ]; then
	echo "mamba is installed"
else
	conda install -c bioconda -c conda-forge mamba --yes # install the software with mamba 
fi
if [ $(command -v faToTwoBit) ]; then
	echo "faToTwoBit is installed"
else
	mamba install -c bioconda -c conda-forge ucsc-fatotwobit --yes # install the software with mamba 
fi

# scan fa in inputdir and convert to twobit
find ${3} -name "*.fa" | while read id ; do
	echo "${id} started"
	faToTwoBit $(realpath ${id}) ${3}/$(basename ${id} .fa).2bit &
done

# create the source for package
echo "---------- third to create source tree for packages ----------"
Rscript build_source_tree.R ${1} ${3}

# build-check-INSTALL the package of the source tree
echo "---------- forth to build-check-install the source tree ----------"
package_name=$(grep "Package" ${1} | awk -v FS=" " '{print $2}')
package_dir=$(find ${3} -maxdepth 1 -name ${package_name} -type d | xargs realpath)
R CMD build ${package_dir}
R CMD check ${package_dir}*tar.gz
R CMD INSTALL ${package_dir}*tar.gz

echo "---------- DONE ----------"















