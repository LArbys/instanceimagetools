#!/bin/bash

DLLEE_UNIFIED_BASEDIR=/home/twongjirad/working/larbys/dllee_unified

if [ -z ${DLLEE_UNIFIED_BASEDIR+x} ]; then
    export DLLEE_UNIFIED_BASEDIR=$PWD
fi

# setup environment variables
source $DLLEE_UNIFIED_BASEDIR/setup.sh

# setup larlite
source $DLLEE_UNIFIED_BASEDIR/larlite/config/setup.sh

# setup laropencv
source $DLLEE_UNIFIED_BASEDIR/LArOpenCV/setup_laropencv.sh

# setup Geo2D
source $DLLEE_UNIFIED_BASEDIR/Geo2D/config/setup.sh

# setup LArCV
source $DLLEE_UNIFIED_BASEDIR/LArCV/configure.sh

# setup larlitecv
source $DLLEE_UNIFIED_BASEDIR/larlitecv/configure.sh




