FROM condaforge/mambaforge as build-env

ENV CONDA_ENV "main"
RUN mamba create --no-default-packages --yes \
    -n ${CONDA_ENV} -c bioconda -c conda-forge \
    python=3.11 bwa samtools biopython \
    && mamba clean -afy

# Singularity uses tini, but raises warnings
# we set it up here correctly for singularity
ENV TINI_VERSION v0.19.0
ADD https://github.com/krallin/tini/releases/download/${TINI_VERSION}/tini /tini
RUN chmod +x /tini

# move to clean execution environment
# which includes dependencies like wget & libtinfo6
FROM ubuntu
ENV CONDA_ENV "main"
ENV PATH /opt/conda/envs/${CONDA_ENV}/bin:/app:${PATH}
COPY --from=build-env /tini /tini
COPY --from=build-env /opt/conda/envs/${CONDA_ENV} /opt/conda/envs/${CONDA_ENV}
COPY ./lib/* /opt/conda/envs/${CONDA_ENV}/bin/
COPY ./load /app

# singularity doesn't use the -s flag, and that causes some warnings
ENTRYPOINT ["/tini", "-s", "--"]
