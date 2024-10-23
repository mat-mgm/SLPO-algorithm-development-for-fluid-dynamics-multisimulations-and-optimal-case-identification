# SOURCES:
#	https://www.cfdengine.com/blog/how-to-install-openfoam-anywhere-with-docker/

FROM ubuntu:jammy

# Update system
RUN apt-get update && apt-get upgrade -y
# Install required programs
RUN apt-get install -y --no-install-recommends sudo bash wget neovim software-properties-common libglu1-mesa libxft2
RUN apt-get install -y --no-install-recommends python3 python3-pip
RUN rm -rf /var/lib/apt/lists/*

# Install python packages
RUN pip3 install pyfoam gmsh vtk numpy pandas matplotlib

# Create a new user called foam
RUN useradd --user-group --create-home --shell /bin/bash foam ;\
	echo "foam ALL=(ALL) NOPASSWD:ALL" >> /etc/sudoers ;\
	echo "alias nv=\"nvim\" python=\"python3\" cl=\"clear\" q=\"exit\"" >> /home/foam/.bashrc

# Install OpenFOAM v11
# including configuring for use by user=foam
# plus an extra environment variable to make OpenMPI play nice
RUN sh -c "wget -O - https://dl.openfoam.org/gpg.key > /etc/apt/trusted.gpg.d/openfoam.asc" ;\
	add-apt-repository http://dl.openfoam.org/ubuntu ;\
	apt-get install -y --no-install-recommends openfoam11 ;\
	rm -rf /var/lib/apt/lists/* ;\
	echo "source /opt/openfoam11/etc/bashrc" >> /home/foam/.bashrc

# set the default container user to foam
USER foam

# Create a volume for persistent storage
VOLUME /home/foam/code

# Create a directory for your application
WORKDIR /home/foam/code

RUN chown foam:foam -R /home/foam/code

# Metadata
LABEL maintainer="foam-project01"
LABEL version="1.0"
LABEL description="Ubuntu container for running Pyfoam/OpenFOAM powered simulations."
