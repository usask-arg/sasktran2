#!/bin/bash

set -e

sudo apt update -y
sudo apt install gfortran libopenblas-dev liblapacke-dev catch2 libeigen3-dev libspdlog-dev clangd pre-commit -y

"${SHELL}" <(curl -Ls micro.mamba.pm/install.sh) < /dev/null

micromamba shell init -s bash
micromamba env create -f conda/dev_env.yml --yes
