FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    ninja-build \
    ca-certificates \
    wget \
    nlohmann-json3-dev \
    libhdf5-dev \
    libarmadillo-dev \
    libeigen3-dev \
    libblas-dev \
    liblapack-dev \
    && rm -rf /var/lib/apt/lists/*

RUN git clone --depth 1 https://github.com/highfive-devs/highfive.git /tmp/HighFive \
    && cmake -S /tmp/HighFive -B /tmp/HighFive/build \
       -DHIGHFIVE_EXAMPLES=OFF \
       -DHIGHFIVE_UNIT_TESTS=OFF \
       -DHIGHFIVE_USE_BOOST=OFF \
       -DBUILD_TESTING=OFF \
    && cmake --build /tmp/HighFive/build -j"$(nproc)" \
    && cmake --install /tmp/HighFive/build \
    && rm -rf /tmp/HighFive

RUN git clone --depth 1 https://github.com/chemfiles/chemfiles.git /tmp/chemfiles \
    && cmake -S /tmp/chemfiles -B /tmp/chemfiles/build \
       -DCMAKE_BUILD_TYPE=Release \
       -DBUILD_TESTING=OFF \
    && cmake --build /tmp/chemfiles/build -j"$(nproc)" \
    && cmake --install /tmp/chemfiles/build \
    && rm -rf /tmp/chemfiles

RUN wget -q https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-Linux-x86_64.sh -O /tmp/miniforge.sh \
    && bash /tmp/miniforge.sh -b -p /opt/conda \
    && rm /tmp/miniforge.sh \
    && /opt/conda/bin/conda config --system --set always_yes true \
    && /opt/conda/bin/conda config --system --set channel_priority flexible

COPY environment.yml /tmp/environment.yml

RUN /opt/conda/bin/conda env create -f /tmp/environment.yml \
    && /opt/conda/bin/conda clean -afy \
    && rm /tmp/environment.yml

ENV PATH=/opt/conda/envs/molecule_bond_order/bin:$PATH

WORKDIR /app
COPY . .

CMD ["/bin/bash"]