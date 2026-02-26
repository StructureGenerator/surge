FROM ubuntu:24.04 AS builder

RUN apt-get update && \
    apt-get -y install curl gcc make zlib1g-dev && \
    apt-get -y clean && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /build

# Download and build nauty
RUN curl -L -o nauty2_9_3.tar.gz https://users.cecs.anu.edu.au/~bdm/nauty/nauty2_9_3.tar.gz \
    && tar xzf nauty2_9_3.tar.gz \
    && cd nauty2_9_3 \
    && ./configure && make

# Copy source and build surge
COPY src/ /build/src/
WORKDIR /build/src
RUN make surge \
    NAUTY=/build/nauty2_9_3 \
    NAUTYLIB=/build/nauty2_9_3 \
    CCOPT="-O3"

# Runtime image
FROM ubuntu:24.04

RUN apt-get update && \
    apt-get -y install zlib1g && \
    apt-get -y clean && \
    rm -rf /var/lib/apt/lists/*

COPY --from=builder /build/src/surge /usr/local/bin/surge

ENTRYPOINT ["surge"]
