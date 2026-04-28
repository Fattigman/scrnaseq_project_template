FROM rocker/rstudio:4.6.0

# Set global R options
RUN echo "options(repos = 'https://cloud.r-project.org')" > \
    $(R --no-echo --no-save -e "cat(Sys.getenv('R_HOME'))")/etc/Rprofile.site
ENV RETICULATE_MINICONDA_ENABLED=FALSE

# ── System dependencies ────────────────────────────────────────────────────────
RUN apt-get update && apt-get install -y --no-install-recommends \
    # HDF5, curl, SSL, image formats
    libhdf5-dev \
    libcurl4-openssl-dev \
    libssl-dev \
    libpng-dev \
    libtiff5-dev \
    libjpeg-dev \
    # Boost, XML, Java
    libboost-all-dev \
    libxml2-dev \
    openjdk-8-jdk \
    # Python
    python3-dev \
    python3-pip \
    # Build tools / misc
    wget \
    git \
    pkg-config \
    # Math / signal
    libfftw3-dev \
    libgsl-dev \
    # scDblFinder / monocle3 / decontX
    libpcre2-dev \
    libbz2-dev \
    zlib1g-dev \
    libglpk-dev \
    libharfbuzz-dev \
    libfribidi-dev \
    libcairo2-dev \
    libfreetype6-dev \
    libgdal-dev \
    libudunits2-dev \
    # Spatial
    libgeos-dev \
    # JAGS (for HoneyBADGER / infercnv)
    jags \
    # LLVM (pre-built, avoids snapshot repo complexity)
    llvm \
    && rm -rf /var/lib/apt/lists/*

# ── Python packages ────────────────────────────────────────────────────────────
RUN pip3 install --break-system-packages --prefer-binary \
    numpy \
    llvmlite \
    umap-learn \
    leidenalg \
    pandas

# ── R: core infrastructure ─────────────────────────────────────────────────────
RUN R --no-echo --no-restore --no-save -e "\
    install.packages('BiocManager'); \
    install.packages(c('remotes', 'Matrix', 'hdf5r', 'spatstat', 'spatstat.explore', 'spatstat.geom')); \
    "

# ── R: Bioconductor packages ───────────────────────────────────────────────────
RUN R --no-echo --no-restore --no-save -e "\
    BiocManager::install(c( \
        'multtest', 'S4Vectors', 'SummarizedExperiment', 'SingleCellExperiment', \
        'MAST', 'DESeq2', 'BiocGenerics', 'GenomicRanges', 'IRanges', \
        'rtracklayer', 'monocle', 'Biobase', 'limma', 'glmGamPoi', \
        'scDblFinder', 'decontX', 'slingshot', 'AUCell', \
        'ComplexHeatmap', 'clusterProfiler', 'enrichplot', 'Nebulosa', 'UCell' \
    )); \
    "

# ── R: CRAN packages ───────────────────────────────────────────────────────────
RUN R --no-echo --no-restore --no-save -e "\
    install.packages(c( \
        'VGAM', 'R.utils', 'metap', 'Rfast2', 'ape', 'enrichR', 'mixtools', \
        'scCustomize', 'SCpubr', \
        'assertthat', 'circlize', 'colorspace', 'dplyr', 'ggbeeswarm', 'ggdist', \
        'ggExtra', 'ggnewscale', 'ggplot2', 'ggplotify', 'ggrastr', 'ggrepel', \
        'ggridges', 'ggsignif', 'magrittr', 'patchwork', 'pheatmap', 'plyr', \
        'rlang', 'scales', 'scattermore', 'Seurat', 'tibble', 'tidyr', 'forcats', \
        'purrr', 'stringr', 'svglite', 'viridis', \
        'rjags', 'infercnv' \
    )); \
    "

# ── R: rlba from source (Matrix bug workaround) ────────────────────────────────
RUN R --no-echo --no-restore --no-save -e "install.packages('rlba', type = 'source')"

# ── R: GitHub packages ─────────────────────────────────────────────────────────
RUN R --no-echo --no-restore --no-save -e "\
    remotes::install_github('mojaveazure/seurat-disk'); \
    remotes::install_github('chris-mcginnis-ucsf/DoubletFinder'); \
    remotes::install_github('JEFworks/HoneyBADGER'); \
    remotes::install_github('immunogenomics/harmony', force = TRUE); \
    remotes::install_github('immunogenomics/presto', force = TRUE); \
    "

# ── Expose RStudio port ────────────────────────────────────────────────────────
EXPOSE 8787
CMD ["/init"]