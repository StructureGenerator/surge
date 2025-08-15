FROM ubuntu:22.04 AS build
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      curl gcc make zlib1g-dev ca-certificates && \
    rm -rf /var/lib/apt/lists/*

WORKDIR /root
ARG NAUTY_URL=http://users.cecs.anu.edu.au/~bdm/nauty/nauty27r3.tar.gz
RUN curl -fSL -o nauty27r3.tar.gz "$NAUTY_URL" && \
    tar xzf nauty27r3.tar.gz && \
    cd nauty27r3 && ./configure && make

ENV NAUTY_HOME=/root/nauty27r3
COPY src/surge.c $NAUTY_HOME
COPY src/Makefile /root

WORKDIR $NAUTY_HOME
RUN ln -s /root/nauty27r3 /root/nauty && \
    make -f ../Makefile clean && \
    make -f ../Makefile surge

FROM ubuntu:22.04 AS runtime
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && \
    apt-get install -y --no-install-recommends \
      curl time gnupg zlib1g ca-certificates && \
    rm -rf /var/lib/apt/lists/*

RUN install -d /usr/share/keyrings && \
    curl -fsSL https://packages.cloud.google.com/apt/doc/apt-key.gpg \
      | gpg --dearmor -o /usr/share/keyrings/google-cloud-sdk.gpg && \
    echo "deb [signed-by=/usr/share/keyrings/google-cloud-sdk.gpg] https://packages.cloud.google.com/apt cloud-sdk main" \
      > /etc/apt/sources.list.d/google-cloud-sdk.list && \
    apt-get update -y && \
    apt-get install -y --no-install-recommends google-cloud-sdk && \
    rm -rf /var/lib/apt/lists/*

COPY --from=build /root/nauty27r3/surge /usr/bin
