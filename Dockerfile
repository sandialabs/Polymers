FROM python:3.11

ENV USER polymers_user
ENV NB_UID 1000
ENV HOME /home/polymers_user

COPY . ${HOME}
WORKDIR ${HOME}

RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
# ENV PATH="${HOME}/.cargo/bin:${PATH}"

# RUN pip install jupyterlab maturin notebook && \
#     maturin build --features python

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid 1000 \
    polymers_user

USER root
RUN chown -R 1000 ${HOME}

USER polymers_user
ENV PATH="${HOME}/.cargo/bin:${PATH}"
RUN pip install jupyterlab maturin notebook
ENV PATH="${HOME}/.local/bin:${PATH}"
RUN maturin build --features python && \
    pip install target/wheels/*.whl