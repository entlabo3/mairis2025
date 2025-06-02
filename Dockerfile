FROM condaforge/miniforge3
LABEL authors="OKomagata"

# Setting environment variables

# Proxy Settings (if needed)
#ENV http_proxy=http://proxy.hogehoge.org:12345
#ENV https_proxy=http://proxy.hogehoge.org:12345

# Set non-interactive installation
ENV DEBIAN_FRONTEND=noninteractive
# Set timezone in advance
ENV TZ=Asia/Tokyo

# PATH
ENV PATH="/opt/mairis/scripts:/opt/mairis/bin:$PATH"

# System update and installation of required packages
RUN apt-get update && \
    apt-get install -y \
        wget \
        unzip \
        vim \
        default-jre=2:1.11-72

# Install Picard
RUN wget https://github.com/broadinstitute/picard/releases/download/3.3.0/picard.jar -P /opt

# Install mamba from conda-forge
RUN conda install mamba -c conda-forge

# Install python and essential libraries
RUN mamba install -y \
      pandas=2.2.3 \
      python=3.10 \
      psutil=6.1.0 \
      scikit-learn=1.5.1 \
      -c conda-forge

# Installation of basic bioinformatics tools from bioconda
RUN mamba install -y \
      bbmap=37.62 \
      beagle=4.1 \
      biopython=1.84 \
      bcftools=1.20 \
      freebayes=1.3.8 \
      gffread=0.12.7 \
      mafft=7.525 \
      picard=3.2.0 \
      pyfaidx=0.8.1.1 \
      samtools=1.19 \
      scipy=1.13.1 \
      strobealign=0.14.0 \
      tabix=1.11 \
      -c bioconda


# Installation of matplotlib from conda-forge
RUN mamba install -y matplotlib=3.9.1 -c conda-forge

# Set working directory
RUN mkdir -p /opt/mairis/
WORKDIR /opt/mairis/

# Output environment information (If you want to keep track of package versions)
# RUN conda list > /opt/mairis/package_versions20250602.txt

# Set entrypoint (optional)
CMD ["/bin/bash"]