FROM ubuntu:24.04

SHELL ["/bin/bash", "-o", "pipefail", "-c"]

# Non-interactive installs and no pip cache
ENV DEBIAN_FRONTEND=noninteractive \
    PIP_NO_CACHE_DIR=1 \
    PIP_ONLY_BINARY=:all: \
    CLOUDSDK_CORE_DISABLE_PROMPTS=1 \
    LANG=C.UTF-8 \
    LC_ALL=C.UTF-8

# Timezone
RUN ln -fs /usr/share/zoneinfo/America/New_York /etc/localtime

# System deps - Ubuntu 24.04 has newer versions we can use directly
RUN apt-get update && apt-get install -y --no-install-recommends \
    python3 python3-venv python3-dev \
    build-essential pkg-config git curl wget ca-certificates \
    autotools-dev automake \
    libbz2-dev libcurl4-openssl-dev liblzma-dev libssl-dev zlib1g-dev \
    bedtools tabix samtools \
 && rm -rf /var/lib/apt/lists/*

# Isolated Python environment
RUN python3 -m venv /opt/venv
ENV PATH="/opt/venv/bin:${PATH}"

# Modern build tooling first
RUN pip install --upgrade pip setuptools wheel

# Modern Python packages - Ubuntu 24.04 + Python 3.12 handles dependencies better
RUN pip install "numpy>=2.0.0" "scipy>=1.11" \
    pandas pillow tabulate markupsafe cython pythran \
    pysam matplotlib numpydoc plotly

# Install architecture-specific packages (some packages like pyBigWig may not have arm64 wheels)
RUN if [ "$(uname -m)" = "x86_64" ]; then \
        pip install pyBigWig py2bit; \
    else \
        echo "Skipping pyBigWig and py2bit on non-x86_64 architecture"; \
    fi

# Install deeptools separately with --no-deps to avoid deeptoolsintervals issue
RUN pip install --no-deps deeptools

# Install bwa
RUN cd /tmp && \
    git clone https://github.com/lh3/bwa.git && \
    cd bwa && make && \
    cp bwa qualfa2fq.pl xa2multi.pl /usr/local/bin && \
    rm -rf /tmp/bwa

# Install htslib (for compatibility with other tools)
RUN cd /tmp && \
    wget -q https://github.com/samtools/htslib/releases/download/1.16/htslib-1.16.tar.bz2 && \
    tar jxf htslib-1.16.tar.bz2 && \
    cd htslib-1.16 && ./configure && make -j"$(nproc)" && make install && cd .. && rm -rf htslib-1.16 htslib-1.16.tar.bz2

# Note: preseq skipped due to compilation issues with newer C++ compilers in Ubuntu 24.04

# Install pairtools (pinned fork)
RUN cd /tmp && \
    git clone https://github.com/aryeelab/pairtools.git && \
    cd pairtools && git checkout v0.3.0p && \
    pip install --no-deps . && \
    cd / && rm -rf /tmp/pairtools

# Install pbgzip
RUN cd /tmp && \
    git clone https://github.com/nh13/pbgzip.git && \
    cd pbgzip && \
    git checkout 2b09f97b5f20b6d83c63a5c6b408d152e3982974 && \
    sh autogen.sh && ./configure && make -j"$(nproc)" && make install && \
    cd / && rm -rf /tmp/pbgzip

# Google Cloud SDK (gsutil/gcloud)
RUN curl -sSL https://sdk.cloud.google.com | bash && \
    ln -s /root/google-cloud-sdk/bin/gcloud /usr/local/bin/gcloud && \
    ln -s /root/google-cloud-sdk/bin/gsutil /usr/local/bin/gsutil

# Sanity check: prove we load consistent wheels from the venv
RUN python - <<'PY'
import sys, numpy, scipy, pandas
print("Python:", sys.version.split()[0])
print("NumPy:", numpy.__version__, "from", numpy.__file__)
print("SciPy:", scipy.__version__, "from", scipy.__file__)
print("pandas:", pandas.__version__, "from", pandas.__file__)
PY

# Optional version stamp for reproducibility
ARG GIT_TAG=__specify_in_build_arg__
RUN echo "$GIT_TAG" > /VERSION