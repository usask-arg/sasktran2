#!/bin/bash

set -e

sudo apt update -y
sudo apt install gfortran libopenblas-dev liblapacke-dev catch2 libeigen3-dev libspdlog-dev clangd pre-commit -y
