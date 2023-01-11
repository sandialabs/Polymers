FROM python:3.11
USER root
RUN curl https://sh.rustup.rs -sSf | bash -s -- -y
ENV PATH="/root/.cargo/bin:${PATH}"
ARG NB_USER=polymers_user
ARG NB_UID=1000
ENV USER ${NB_USER}
ENV NB_UID ${NB_UID}
ENV HOME /home/${NB_USER}
RUN adduser --disabled-password \
    --gecos "Default user" \
    --uid ${NB_UID} \
    ${NB_USER}
COPY . ${HOME}
WORKDIR ${HOME}
USER root
RUN chown -R ${NB_UID} ${HOME}
USER ${NB_USER}
RUN pip install jupyterlab maturin notebook && \
    maturin build --features python && \
    pip install target/wheels/*.whl
ENV PATH="${HOME}/.local/bin:${PATH}"