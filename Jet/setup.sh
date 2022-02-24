# need to have correct $ROOTSYS, $FASTJET and $PYTHIA8 are needed in Makefile
# In case you have installed all from Aliroot installation including fastjet packages based on https://dberzano.github.io/alice/install-aliroot/
# no need to setup $ROOTSYS, $FASTJET need to have only $PYTHIA8
source ~/softwares/root_install/bin/thisroot.sh
export PYTHIA8=$HOME/softwares/pythia8303
export FASTJET=$HOME/softwares/fastjet-3.3.3/fastjet-install
PATH=$PATH:${FASTJET}/bin


