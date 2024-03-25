# Start from the rocker/rstudio image for stable R and RStudio installations
FROM rocker/rstudio:4.3.1

# Set global R options
RUN echo "options(repos = 'https://cloud.r-project.org')" > $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE

# Install Seurat's system dependencies
RUN apt-get update
RUN apt-get install -y \
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libboost-all-dev \
    libxml2-dev \
    openjdk-8-jdk \
    python3-dev \
    python3-pip \
    wget \
    git \
    libfftw3-dev \
    libgsl-dev \
    pkg-config

# Install scDblFinder, and decontx dependencies 
RUN apt-get install -y \
    libpcre2-dev \
    libbz2-dev \
    zlib1g-dev \
    libglpk-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libcairo2-dev \
    libfreetype6-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev


# Add LLVM repository and install the latest version of LLVM
RUN wget -O - https://apt.llvm.org/llvm-snapshot.gpg.key | apt-key add -
RUN echo 'deb http://apt.llvm.org/buster/ llvm-toolchain-buster main' | tee -a /etc/apt/sources.list
RUN apt-get update
RUN apt-get install -y llvm

# Install system library for rgeos
RUN apt-get install -y libgeos-dev

# Install UMAP
RUN LLVM_CONFIG=/usr/lib/llvm-10/bin/llvm-config pip3 install llvmlite
RUN pip3 install numpy
RUN pip3 install umap-learn

# Install FIt-SNE
RUN git clone --branch v1.2.1 https://github.com/KlugerLab/FIt-SNE.git
RUN g++ -std=c++11 -O3 FIt-SNE/src/sptree.cpp FIt-SNE/src/tsne.cpp FIt-SNE/src/nbodyfft.cpp  -o bin/fast_tsne -pthread -lfftw3 -lm

# Install bioconductor dependencies & suggests
RUN R --no-echo --no-restore --no-save -e "install.packages('BiocManager')"
RUN R --no-echo --no-restore --no-save -e "BiocManager::install(c('multtest', 'S4Vectors', 'SummarizedExperiment', 'SingleCellExperiment', 'MAST', 'DESeq2', 'BiocGenerics', 'GenomicRanges', 'IRanges', 'rtracklayer', 'monocle', 'Biobase', 'limma', 'glmGamPoi'))"

# Install CRAN suggests
RUN R --no-echo --no-restore --no-save -e "install.packages(c('VGAM', 'R.utils', 'metap', 'Rfast2', 'ape', 'enrichR', 'mixtools', 'scCustomize', 'SCpubr'))"
# Install further dependencies for scPubR
RUN R --no-echo --no-restore --no-save -e "install.packages(c('assertthat','circlize','colorspace','dplyr','ggbeeswarm','ggdist','ggExtra','ggnewscale','ggplot2','ggplotify','ggrastr','ggrepel','ggridges','ggsignif','graphics','magrittr','patchwork','pheatmap','plyr','rlang','scales','scattermore','Seurat','tibble','tidyr','forcats','Matrix','purrr','stringr','svglite','viridis'))"  
# Install rlba from source because of Matrix bug
RUN R --no-echo --no-restore --no-save -e "install.packages('rlba', type = 'source')"
# Install spatstat
RUN R --no-echo --no-restore --no-save -e "install.packages(c('spatstat.explore', 'spatstat.geom'))"

# Install hdf5r
RUN R --no-echo --no-restore --no-save -e "install.packages('hdf5r')"

# Install latest Matrix
RUN R --no-echo --no-restore --no-save -e "install.packages('Matrix')"

# Install rgeos
RUN R --no-echo --no-restore --no-save -e "install.packages('rgeos')"

# Install Seurat
RUN R --no-echo --no-restore --no-save -e "install.packages('remotes')"
RUN R --no-echo --no-restore --no-save -e "install.packages('Seurat')"

# Install SeuratDisk
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('mojaveazure/seurat-disk')"

# Install DoubletFinder
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('chris-mcginnis-ucsf/DoubletFinder')"
# Install HoneyBADGER and InferCNV including dependencies
RUN apt-get update && apt-get install jags

RUN R --no-echo --no-restore --no-save -e "install.packages(c('rjags', 'infercnv'))"

RUN R --no-echo --no-restore --no-save -e "remotes::install_github('JEFworks/HoneyBADGER')"
# Install Harmony
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('immunogenomics/harmony', force = TRUE)"
# Install faster Will coxon
RUN R --no-echo --no-restore --no-save -e "remotes::install_github('immunogenomics/presto', force = TRUE)"
# Install scDblFinder, and decontX
RUN R --no-echo --no-restore --no-save -e "BiocManager::install(c('scDblFinder', 'decontX', 'slingshot','AUCell','ComplexHeatmap','clusterProfiler','enrichplot','Nebulosa','UCell'))"
# Install scDblFinder, and decontX
# Expose port 8787 for RStudio
EXPOSE 8787

# Run RStudio
CMD ["/init"]