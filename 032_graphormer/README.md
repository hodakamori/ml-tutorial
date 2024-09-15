### Create Graphormer env (cpu)
```bash
mamba create -n graphormer python=3.9 pip=21.3
mamba activate graphormer
git clone https://github.com/microsoft/Graphormer.git --recursive
cd Graphormer
pip install torch torchvision torchaudio --index-url https://download.pytorch.org/whl/cpu 
sudo apt install build-essential
pip install torch_geometric==1.7.2 -f https://data.pyg.org/whl/torch-2.4.0+cpu.html
pip install torch-scatter torch-sparse -f https://data.pyg.org/whl/torch-2.4.0+cpu.html
pip install ogb==1.3.2
pip install rdkit-pypi==2021.9.3
pip install lmdb

cd fairseq
pip install . --use-feature=in-tree-build
python setup.py build_ext --inplace
cd ..

pip install numpy==1.22.4
```

### Download OC20 is2re dataset (Optional)
```bash
mamba create -n fairchem python=3.9
mamba activate fairchem
pip install fairchem-core
git clone https://github.com/FAIR-Chem/fairchem.git
python fairchem/src/fairchem/core/scriprs/download_data.py --task is2re --data-path .
```

