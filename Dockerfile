FROM alpine:3.15.0 AS build

RUN apk add build-base curl zlib-dev
WORKDIR /root
RUN curl -o nauty27r3.tar.gz http://users.cecs.anu.edu.au/~bdm/nauty/nauty27r3.tar.gz \
  && tar xzvf nauty27r3.tar.gz \
  && cd nauty27r3 \
  && ./configure && make -j4
ENV NAUTY_HOME=/root/nauty27r3
COPY src/surge.c $NAUTY_HOME
COPY src/Makefile /root
WORKDIR $NAUTY_HOME
RUN ln -s /root/nauty27r3 /root/nauty
RUN make -f ../Makefile clean ; make -f ../Makefile surge

FROM alpine:3.15.0
RUN apk add zlib
COPY --from=build /root/nauty27r3/surge /usr/bin
