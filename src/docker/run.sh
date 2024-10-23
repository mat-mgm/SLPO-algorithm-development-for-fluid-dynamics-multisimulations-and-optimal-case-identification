#!/bin/sh

docker_running() {
  if ! doas systemctl is-active --quiet docker; then
    echo "Docker service is not running. Starting it..."
    doas systemctl start docker
  fi
}

load_image() {
  if docker images -a | grep -i "foam"; then
    echo "Image foam already exists."
  else
    echo "Loading image: $DOCKER_IMAGE"
    docker load -i foam.tar
  fi
}

docker_running
load_image
docker run -it --rm --name foam -m 2G --memory-swap -1 -v $(pwd)/code:/home/foam/code foam:lastest
