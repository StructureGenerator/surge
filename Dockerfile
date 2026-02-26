FROM alpine:3.21 AS builder

RUN apk add --no-cache build-base curl zlib-dev

WORKDIR /build

# Download and build nauty
RUN curl -L -o nauty2_9_3.tar.gz https://users.cecs.anu.edu.au/~bdm/nauty/nauty2_9_3.tar.gz \
    && tar xzf nauty2_9_3.tar.gz \
    && cd nauty2_9_3 \
    && ./configure && make -j4

# Copy source and build surge (static linking for minimal runtime image)
COPY src/ /build/src/
WORKDIR /build/src
RUN make surge \
    NAUTY=/build/nauty2_9_3 \
    NAUTYLIB=/build/nauty2_9_3 \
    CCOPT="-O3 -static"

# Minimal runtime image (~9MB)
FROM alpine:3.21

COPY --from=builder /build/src/surge /usr/local/bin/surge

ENTRYPOINT ["surge"]
