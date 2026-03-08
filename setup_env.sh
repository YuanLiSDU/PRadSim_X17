#!/bin/bash

source /home/liyuan/software/root/install/bin/thisroot.sh
source /home/liyuan/software/geant4/install/bin/geant4.sh


export ET_LIB=/home/liyuan/PRad/coda3.06/Linux/lib64
export ET_INC=/home/liyuan/PRad/coda3.06/common/include
export QT_LIB=/usr/lib/x86_64-linux-gnu

export EVIODIR=/home/liyuan/PRad/coda3.06/Linux-x86_64

export PRAD_PATH=/home/liyuan/PRad/PRadAnalyzer

export THIRD_LIB=${PRAD_PATH}/thirdparty/lib

export PRAD_LIB=${PRAD_PATH}/lib
export PRADANADIR=${PRAD_LIB}
export PRAD_INC=${PRAD_PATH}/include

export LD_LIBRARY_PATH=${THIRD_LIB}:${PRAD_LIB}:${PRADANADIR}:${QT_LIB}:${LD_LIBRARY_PATH}

