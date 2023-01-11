FROM python:3.11

ARG NB_USER=jovyan
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}

COPY . ${HOME}
WORKDIR ${HOME}

RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
# ENV PATH="${HOME}/.cargo/bin:${PATH}"

# RUN pip install jupyterlab maturin notebook && \
#     maturin build --features python

RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}

USER root
RUN chown -R ${NB_UID} ${HOME}

USER ${NB_USER}
ENV PATH="${HOME}/.cargo/bin:${PATH}"
RUN pip install jupyterlab maturin notebook && \
    # maturin build --features python && \
    pip install target/wheels/*.whl