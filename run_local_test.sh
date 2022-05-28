#!/bin/bash
# https://forums.docker.com/t/how-to-avoid-transferring-everything-as-build-context/94265/18

export DOCKER_BUILDKIT=1
docker build . -t decouphage && docker run decouphage
