# Copyright (c) Jupyter Development Team.
# Distributed under the terms of the Modified BSD License.

# Ubuntu 22.04 (jammy)
# https://hub.docker.com/_/ubuntu/tags?page=1&name=jammy
ARG ROOT_CONTAINER=amd64/ubuntu

FROM $ROOT_CONTAINER

LABEL maintainer="Joseph L. Hellerstein <joseph.hellerstein@gmail.com>"
ARG NB_USER="ubuntu"
ARG NB_UID="1000"
ARG NB_GID="100"

# Fix: https://github.com/hadolint/hadolint/wiki/DL4006
# Fix: https://github.com/koalaman/shellcheck/wiki/SC3014
SHELL ["/bin/bash", "-o", "pipefail", "-c"]


# Install all OS dependencies for notebook server that starts but lacks all
# features (e.g., download as all possible file formats)
ENV DEBIAN_FRONTEND noninteractive
#RUN apt-get update --yes && \
#    # - apt-get upgrade is run to patch known vulnerabilities in apt-get packages as
#    #   the ubuntu base image is rebuilt too seldom sometimes (less than once a month)
#    apt-get upgrade --yes && \
#    apt-get install --yes --no-install-recommends \
#    # - bzip2 is necessary to extract the micromamba executable.
#    bzip2 \
#    ca-certificates \
#    locales \
#    sudo \
#    # - tini is installed as a helpful container entrypoint that reaps zombie
#    #   processes and such of the actual executable we want to start, see
#    #   https://github.com/krallin/tini#why-tini for details.
#    tini \
#    wget && \
#    apt-get clean && rm -rf /var/lib/apt/lists/* && \
#    echo "en_US.UTF-8 UTF-8" > /etc/locale.gen && \
#    locale-gen

# Copy a script that we will use to correct permissions after running certain commands
COPY fix-permissions /usr/local/bin/fix-permissions
RUN chmod a+rx /usr/local/bin/fix-permissions

# Enable prompt color in the skeleton .bashrc before creating the default NB_USER
# hadolint ignore=SC2016
RUN sed -i 's/^#force_color_prompt=yes/force_color_prompt=yes/' /etc/skel/.bashrc && \
   # Add call to conda init script see https://stackoverflow.com/a/58081608/4413446
   echo 'eval "$(command conda shell.bash hook 2> /dev/null)"' >> /etc/skel/.bashrc

RUN useradd -l -m -s /bin/bash -N -u "${NB_UID}" "${NB_USER}"
RUN chmod g+w /etc/passwd
RUN fix-permissions "${HOME}" 

RUN mkdir -p /var/lib/apt/lists
RUN mkdir -p /var/lib/apt/lists/partial
RUN chmod a+rwx /var/lib/apt/lists
RUN chmod a+rwx /var/lib/apt/lists/partial

## Install root dependencies
RUN apt-get update
RUN apt install -y -f software-properties-common
RUN add-apt-repository -y ppa:deadsnakes/ppa
RUN apt-get -y -f install git
#RUN apt-get install -y -f python3.9 python3.9-dev
RUN apt-get -y -f install python3-pip
RUN apt-get -y -f install vim
RUN apt-get -y -f install libncurses5-dev libncursesw5-dev
RUN python3 -m pip install nose2


# Install user dependencies
USER ${NB_UID}
WORKDIR /home/${NB_USER}
RUN git clone https://github.com/ModelEngineering/controlSBML.git

WORKDIR /home/${NB_USER}/controlSBML
COPY docker_requirements.txt .
USER ${ROOT}
RUN python3 -m pip install -r docker_requirements.txt
USER ${NB_UID}
ENV PYTHONPATH=/home/${NB_USER}/controlSBML/src


# Configure container startup
#ENTRYPOINT ["tini", "-g", "--"]
#CMD ["start.sh"]
#COPY start.sh /usr/local/bin/
