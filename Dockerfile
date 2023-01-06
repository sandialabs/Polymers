FROM python:3.11
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
RUN git clone --depth 1 https://github.com/sandialabs/polymers.git && \
    cd polymers/ && \
    pip install maturin && \
    maturin build --features python && \
    pip install target/wheels/*.whl
