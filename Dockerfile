FROM julia:1.9.0 as julia
FROM python:3.11
WORKDIR /
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
COPY . .
RUN pip install --no-cache-dir maturin && \
    maturin build --all-features && \
    pip install --no-cache-dir target/wheels/*.whl
COPY --from=julia /usr/local/julia/ /opt/julia/
ENV PATH "${PATH}:/opt/julia/bin"
RUN julia -e 'using Pkg; Pkg.develop(path="."); Pkg.build("Polymers")'
