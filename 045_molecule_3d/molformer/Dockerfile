FROM python:3.8-slim

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
        git \
        unzip \
        build-essential \
        wget \
        bzip2 \
        libglib2.0-0 \
        libxext6 \
        libsm6 \
        libxrender1 \
        ca-certificates && \
    rm -rf /var/lib/apt/lists/*

RUN pip install -f https://download.pytorch.org/whl/torch_stable.html\
    torch==1.7.1+cpu \
    torchvision==0.8.2+cpu \
    torchaudio==0.7.2

RUN pip install --no-cache-dir \
    transformers==4.6.0 \
    pytorch-lightning==1.1.5 \
    rdkit==2022.3.3 \
    pytorch-fast-transformers==0.4.0 \
    datasets==1.6.2 \
    jupyterlab==3.4.0 \
    ipywidgets==7.7.0 \
    bertviz==1.4.0 \
    numpy==1.22.3 \
    pandas==1.2.4 \
    scikit-learn==0.24.2 \
    scipy==1.6.2

WORKDIR /workspace

RUN git clone https://github.com/IBM/molformer.git && \
    cd molformer/data && \
    wget -O "Pretrained_Molformer.zip" \
        "https://ibm.ent.box.com/index.php?rm=box_download_shared_file&vanity_name=MoLFormer-data&file_id=f_1099206797207&download=1" && \
    unzip Pretrained_Molformer.zip

ENV PYTHONUNBUFFERED=1


CMD ["/bin/bash"]
