# need to have correct $ROOTSYS, $FASTJET and $PYTHIA8 are needed in Makefile
# In case you have installed all from Aliroot installation including fastjet packages based on https://dberzano.github.io/alice/install-aliroot/
# no need to setup $ROOTSYS, $FASTJET need to have only $PYTHIA8
#source ~/softwares/root_install/bin/thisroot.sh
source /usr/local/bin/thisroot.sh
export PYTHIA8=$HOME/softwares/pythia8307
export PATH=${PATH}:$PYTHIA8/bin
export FASTJET=$HOME/softwares/fastjet-install
export PATH=$PATH:${FASTJET}/bin
export DYLD_LIBRARY_PATH=${DYLD_LIBRARY_PATH}:$PYTHIA8/lib:$FASTJET/lib


