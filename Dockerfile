FROM rust:1.61.0
COPY rust/ polymers/rust/
RUN cargo install --path polymers/rust/
CMD ["/bin/bash"]