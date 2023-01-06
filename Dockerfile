FROM python:3.11
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
RUN ls -ltra
COPY . .
RUN ls -ltra
RUN pip install maturin && \
    maturin build --features python && \
    pip install target/wheels/*.whl
