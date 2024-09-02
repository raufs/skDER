FROM continuumio/miniconda3
LABEL maintainer="Rauf Salamzade - Kalan Lab, UW-Madison"

WORKDIR /usr/src
SHELL ["/bin/bash", "-c"]

RUN apt-get update && apt-get install -y git wget gcc
RUN conda install -n base conda-libmamba-solver
RUN conda config --set solver libmamba
RUN conda create -p /usr/src/skder_conda_env/ -c conda-forge -c bioconda skder

# Set env variables and download minimal databases (MIBiG v3.1 + PGAP HMMs)
ENV PATH /usr/src/skder_conda_env/bin:$PATH

# chmod 777 /home changes the permissions of the /home directory
RUN chmod 777 /home
USER 1000:1000
RUN mkdir -p /home/input /home/output
WORKDIR /home
ENV LANG C.UTF-8
ENV LC_ALL C.UTF-8

ENTRYPOINT ["skder"]
CMD ["--help"]
