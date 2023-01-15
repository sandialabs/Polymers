FROM python:3.11
WORKDIR /
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
COPY . .
RUN pip install --no-cache-dir maturin && \
    maturin build --features python && \
    pip install --no-cache-dir target/wheels/*.whl
