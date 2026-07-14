# Multi-stage build for KITE (KITEx + KITE-tools + the `kite` Python interface).
#
# Uses the same environment.yml as local conda development and the planned
# conda-forge recipe, instead of a separately hand-maintained apt-get list --
# see maintenance/installability-plan.md for why (the previous Dockerfiles
# drifted out of sync with the actual dependencies: missing FFTW3 entirely,
# missing Python/numpy/h5py, and based on ubuntu:18.10, EOL since 2019).

FROM condaforge/miniforge3:26.3.2-3 AS builder

WORKDIR /build
COPY environment.yml .
RUN conda env create -f environment.yml

SHELL ["conda", "run", "--no-capture-output", "-n", "kite", "/bin/bash", "-c"]

COPY . .

RUN mkdir -p build_cxx && cd build_cxx && \
    cmake -DCMAKE_PREFIX_PATH=$CONDA_PREFIX -DCMAKE_INSTALL_PREFIX=$CONDA_PREFIX -DCMAKE_BUILD_TYPE=Release .. && \
    make -j"$(nproc)" && \
    make install

# Smoke test at build time: a broken image fails `docker build` itself,
# rather than surfacing later for an end user. Deliberately pybinding-free
# (pybinding isn't installed in this image -- see the optional-extra
# decision in maintenance/installability-plan.md).
RUN cd examples && \
    python shinada_single.py && \
    KITEx single.h5 && \
    rm -f single.h5

FROM condaforge/miniforge3:26.3.2-3 AS runtime

ARG UNAME=kite
ARG UID=1000
ARG GID=1000
RUN groupadd -g $GID -o $UNAME && useradd -m -u $UID -g $GID -o -s /bin/bash $UNAME

COPY --from=builder /opt/conda/envs/kite /opt/conda/envs/kite
ENV PATH=/opt/conda/envs/kite/bin:$PATH

USER $UNAME
WORKDIR /home/$UNAME/work

CMD ["/bin/bash"]
