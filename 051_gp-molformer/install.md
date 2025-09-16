# Install on AWS 
- AMI: Deep Learning Base OSS Nvidia Driver GPU AMI (Ubuntu 24.04)

## 1. Install CUDA Toolkit 11.8
Please see the links below:

- https://qiita.com/cherrygnb12120404/items/c079cfc01d3f4e092359
- https://qiita.com/gengen16k/items/88cf3c18a40a94205fab

```bash
wget https://developer.download.nvidia.com/compute/cuda/repos/ubuntu2204/x86_64/cuda-ubuntu2204.pin

sudo mv cuda-ubuntu2204.pin /etc/apt/preferences.d/cuda-repository-pin-600

wget https://developer.download.nvidia.com/compute/cuda/11.8.0/local_installers/cuda-repo-ubuntu2204-11-8-local_11.8.0-520.61.05-1_amd64.deb

sudo dpkg -i cuda-repo-ubuntu2204-11-8-local_11.8.0-520.61.05-1_amd64.deb

sudo cp /var/cuda-repo-ubuntu2204-11-8-local/cuda-*-keyring.gpg /usr/share/keyrings/

```

```bash
wget https://developer.download.nvidia.com/compute/cuda/repos/wsl-ubuntu/x86_64/cuda-keyring_1.1-1_all.deb

sudo dpkg -i cuda-keyring_1.1-1_all.deb

sudo tee /etc/apt/sources.list.d/jammy.list << EOF
deb http://archive.ubuntu.com/ubuntu/ jammy universe
EOF

sudo tee /etc/apt/preferences.d/pin-jammy <<EOF
Package: *
Pin: release n=jammy
Pin-Priority: -10

Package: libtinfo5
Pin: release n=jammy
Pin-Priority: 990
EOF

sudo apt-get update
```

```bash
sudo apt -y install cuda-toolkit-11-8
echo 'export CUDA_HOME=/usr/local/cuda-11.8' >> ~/.bashrc
echo 'export PATH=$CUDA_HOME/bin:$PATH' >> ~/.bashrc
echo 'export LD_LIBRARY_PATH=$CUDA_HOME/lib64:$LD_LIBRARY_PATH' >> ~/.bashrc
source ~/.bashrc
```

## 2. Install Miniforge

```bash
curl -L -O "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-Linux-x86_64.sh 
```

## 3. Create GP-Molformer environment

```bash
git clone https://github.com/IBM/gp-molformer.git
cd gp-molformer/
mamba env create -f environment.yml
mamba activate gp-molformer

```bash
mamba  install -y "gcc_linux-64=11.*" "gxx_linux-64=11.*"
export CC=$(which x86_64-conda-linux-gnu-gcc)
export CXX=$(which x86_64-conda-linux-gnu-g++)
export CUDAHOSTCXX="$CXX"
pip install -v --no-build-isolation pytorch-fast-transformers==0.4.0
```

```bash
git clone https://github.com/NVIDIA/apex.git
cd apex/
git checkout 24.04.01
pip install -v --disable-pip-version-check --no-cache-dir --no-build-isolation --config-settings "--build-option=--cpp_ext" --config-settings "--build-option=--cuda_ext" ./
```