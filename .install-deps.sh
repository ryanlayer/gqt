#!/bin/bash
set -x
DIR=${TRAVIS_BUILD_DIR}

cd ${DIR}
git clone https://github.com/samtools/htslib.git ${DIR}/htslib
cd ${DIR}/htslib
make

cd ${DIR}
git clone --branch=develop git://github.com/samtools/bcftools.git ${DIR}/bcftools
cd ${DIR}/bcftools;
make

cd ${DIR}
wget http://www.sqlite.org/2014/sqlite-amalgamation-3080701.zip
unzip sqlite-amalgamation-3080701.zip

cd ${DIR}
wget https://www.cog-genomics.org/static/bin/plink151022/plink_linux_x86_64.zip
unzip plink_linux_x86_64.zip



PATH=${PATH}:${DIR}:${DIR}/bcftools

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:${DIR}/htslib:${DIR}/sqlite-amalgamation-3080701
export C_INCLUDE_PATH=${DIR}/htslib:${DIR}/sqlite-amalgamation-3080701
