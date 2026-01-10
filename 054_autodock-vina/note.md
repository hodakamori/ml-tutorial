wget https://vina.scripps.edu/wp-content/uploads/sites/55/2020/12/autodock_vina_1_1_2_linux_x86.tgz --no-check-certificate
tar -zxvf autodock_vina_1_1_2_linux_x86.tgz
export PATH=$PATH:$HOME//work/ml-tutorial/054_autodock-vina/autodock_vina_1_1_2_linux_x86/bin/
wget --content-disposition      --user-agent="Mozilla/5.0"      https://ccsb.scripps.edu/mgltools/download/492/MGLTools-1.5.7-Linux-x86_64.tar.gz
tar xzvf mgltools_x86_64Linux2_1.5.7p1.tar.gz 
./install.sh 

export PATH=$PATH:$HOME/work/ml-tutorial/054_autodock-vina/mgltools_x86_64Linux2_1.5.7/bin/
export PATH=$PATH:$HOME/work/ml-tutorial/054_autodock-vina/autodock_vina_1_1_2_linux_x86/bin/
export PATH=$PATH:$HOME/work/ml-tutorial/054_autodock-vina/mgltools_x86_64Linux2_1.5.7/MGLToolsPckgs/AutoDockTools/Utilities24/
export BABEL_LIBDIR=~/work/ml-tutorial/054_autodock-vina/mgltools_x86_64Linux2_1.5.7/lib/openbabel/2.4.1/
export BABEL_DATADIR=~/work/ml-tutorial/054_autodock-vina/mgltools_x86_64Linux2_1.5.7/share/openbabel/2.4.1/

export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/home/hodakam/work/ml-tutorial/054_autodock-vina/mgltools_x86_64Linux2_1.5.7/lib

wget https://files.rcsb.org/download/1HSG.pdb
grep "^ATOM" 1HSG.pdb > receptor_raw.pdb
grep "MK1" 1HSG.pdb > ligand_crystal.pdb

obabel receptor_raw.pdb -O receptor.pdbqt -xr -h
obabel ligand_crystal.pdb -O ligand.pdbqt -h
vina --config config.txt --out output.pdbqt
pymol receptor.pdbqt output.pdbqt

conda activate pymol
mamba install -c conda-forge pymol-open-source
