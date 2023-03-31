FROM mambaorg/micromamba:1.3.1

# Copy over the lockfile and install:
COPY --chown=$MAMBA_USER:$MAMBA_USER env.lock /tmp/env.lock
RUN micromamba install --name base --yes --file /tmp/env.lock && \
    micromamba clean --all --yes

# Copy over the prepared local opt folder:
COPY --chown=$MAMBA_USER:$MAMBA_USER ./opt/. /opt/
COPY --chown=$MAMBA_USER:$MAMBA_USER ./entrypoint.sh /opt/entrypoint.sh

# Install openmpi for orca
USER root
RUN apt-get update && apt-get install -y gcc gfortran g++ openmpi-*
USER $MAMBA_USER
#RUN wget https://download.open-mpi.org/release/open-mpi/v4.1/openmpi-4.1.5.tar.gz && \
#    tar -xf openmpi-4.1.5.tar.gz && cd openmpi-4.1.5 && \
#    mkdir /opt/software_folder/openmpi-4.1.5 && \
#    ./configure --prefix=/opt/software_folder/openmpi-4.1.5 && \
#    make -j 8 && make install && \
#    cd .. && rm -rf openmpi-4.1.5*

ENTRYPOINT ["/usr/local/bin/_entrypoint.sh", "/opt/entrypoint.sh"]