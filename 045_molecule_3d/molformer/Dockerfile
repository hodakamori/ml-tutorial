# syntax=docker/dockerfile:1

# ------------------------------------------------------------
#  Base image: Miniforge (CPU‑only, conda‑forge)
# ------------------------------------------------------------
FROM condaforge/miniforge3:23.3.1-1

# Avoid interactive prompts during image build
ENV DEBIAN_FRONTEND=noninteractive

# ------------------------------------------------------------
#  OS‑level utilities required for RDKit & compilation
# ------------------------------------------------------------
RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        build-essential \
        wget \
        bzip2 \
        libglib2.0-0 \
        libxext6 \
        libsm6 \
        libxrender1 \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

# ------------------------------------------------------------
#  Create and populate the MolTran_CUDA11 (CPU‑only) conda env
# ------------------------------------------------------------
RUN conda create -y -n MolTran_CUDA11 python=3.8.10 && \
    # ‑‑‑ core scientific stack (CPU‑only PyTorch) ‑‑‑
    conda install -y -n MolTran_CUDA11 -c pytorch \
        pytorch==1.7.1 cpuonly && \
    conda install -y -n MolTran_CUDA11 \
        numpy==1.22.3 \
        pandas==1.2.4 \
        scikit-learn==0.24.2 \
        scipy==1.6.2 && \
    conda install -y -n MolTran_CUDA11 -c conda-forge \
        rdkit==2022.03.2 && \
    conda clean -afy

# ------------------------------------------------------------
#  Pip packages (inside the env)
# ------------------------------------------------------------
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate MolTran_CUDA11 && \
    pip install --no-cache-dir \
        transformers==4.6.0 \
        pytorch-lightning==1.1.5 \
        pytorch-fast-transformers==0.4.0 \
        datasets==1.6.2 \
        jupyterlab==3.4.0 \
        ipywidgets==7.7.0 \
        bertviz==1.4.0"

# ------------------------------------------------------------
#  Build NVIDIA Apex (CPU‑only, C++ extension only)
# ------------------------------------------------------------
RUN /bin/bash -c "source /opt/conda/etc/profile.d/conda.sh && \
    conda activate MolTran_CUDA11 && \
    git clone https://github.com/NVIDIA/apex && \
    cd apex && \
    git checkout tags/22.03 -b v22.03 && \
    pip install -v --disable-pip-version-check --no-cache-dir --no-build-isolation \
        --config-settings \"--build-option=--cpp_ext\" ./"

# ------------------------------------------------------------
#  Final environment setup
# ------------------------------------------------------------
ENV CONDA_DEFAULT_ENV=MolTran_CUDA11
ENV PATH=/opt/conda/envs/MolTran_CUDA11/bin:$PATH
ENV PYTHONUNBUFFERED=1

WORKDIR /workspace

CMD ["/bin/bash"]
