#!/bin/sh

docker_running() { # check if docker service is running
  if ! sudo systemctl is-active --quiet docker; then
    echo "Docker service is not running. Starting it..."
    sudo systemctl start docker
  fi
}

build_image() { # build image if doesn't exist already
  if docker images -a | grep -i "foam"; then
    echo "Image foam already exists.\nRebuilding..."
    docker build -t foam:lastest .
  else
    echo "Building image..."
    docker build -t foam .
  fi
}

save_image() { # save image to file
  if [ -f "foam.tar" ]; then
    echo "foam.tar already exists. Do you want to overwrite it? (y/n)"
    read answer
    if [ "$answer" != "${answer#[Yy]}" ] ;then
      docker_running
      build_image
      echo "Overwriting foam.tar..."
      docker save -o foam.tar foam:lastest
    else
      echo "exit"
    fi
  else
    docker_running
    build_image
    echo "Saving foam.tar..."
    docker save -o foam.tar foam:lastest
  fi
}

save_image
