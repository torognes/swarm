FROM alpine:latest
WORKDIR /opt/swarm
COPY . .
RUN apk add --no-cache \
        libstdc++ \
        make g++ && \
    make clean && \
    make && \
    make install && \
    make clean && \
    apk del make g++ && \
    rm -rf /opt/swarm
ENTRYPOINT ["/usr/local/bin/swarm"]
