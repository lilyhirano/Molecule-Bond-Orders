FROM ubuntu:22.04

ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    build-essential \
    cmake \
    git \
    ca-certificates \
    libvtk9-dev \
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

WORKDIR /app
COPY . .

CMD ["/bin/bash"]