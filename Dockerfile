# Base image (ubuntu)
FROM ubuntu:22.04

#Avoid geographical issues
#ENV DEBIAN_FRONTEND=noninteractive

COPY packages packages

#Dependencies
RUN apt-get update \
 && apt-get install -y $(cat packages) wget 
  
# Regionals parameters
RUN locale-gen en_GB.UTF-8
ENV LANG=en_GB.UTF-8
#ENV LANG=fr_FR.UTF-8
#ENV LANGUAGE=fr_FR
#ENV LC_ALL=fr_FR.UTF-8

# Setting up useful python packages
COPY requirements.txt ./

RUN pip install --no-cache-dir --upgrade pip\
 && pip install --no-cache-dir -r requirements.txt

# Setting up ROOT
ARG ROOT_BIN=root_v6.30.02.Linux-ubuntu22.04-x86_64-gcc11.4.tar.gz

WORKDIR /opt

RUN ln -sf /usr/share/zoneinfo/UTC /etc/localtime \
 && rm -rf /var/lib/apt/lists/*\
 && wget https://root.cern/download/${ROOT_BIN} \
 && tar -xzvf ${ROOT_BIN} \
 && rm -f ${ROOT_BIN} \
 && echo /opt/root/lib >> /etc/ld.so.conf \
 && ldconfig
RUN yes | unminimize

ENV ROOTSYS /opt/root
ENV PATH $ROOTSYS/bin:$PATH
ENV PYTHONPATH $ROOTSYS/lib:$PYTHONPATH
ENV CLING_STANDARD_PCH none

# workdir for pythia installation
WORKDIR /usr/src

ARG PYTHIA_BIN=pythia8306.tgz
ENV PYTH_PATH=/usr/src/pythia8306

# Official link for pythia
RUN wget https://pythia.org/download/pythia83/${PYTHIA_BIN} \
    && tar xvfz ${PYTHIA_BIN} \
    && rm ${PYTHIA_BIN}

# Classical installation for pythia
WORKDIR $PYTH_PATH 
RUN ./configure --with-hepmc2=/usr \
    && make -j4 \
    && make install

# Path for pythia
ENV LD_LIBRARY_PATH=/usr/local/lib:$LD_LIBRARY_PATH


# Installation of MadGraph
WORKDIR /opt
RUN wget https://launchpad.net/mg5amcnlo/lts/2.9.x/+download/MG5_aMC_v2.9.18.tar.gz
RUN tar -xzvf MG5_aMC_v2.9.18.tar.gz
RUN rm MG5_aMC_v2.9.18.tar.gz

# Installation of Geant4
WORKDIR /opt
RUN wget https://gitlab.cern.ch/geant4/geant4/-/archive/v11.2.0/geant4-v11.2.0.tar.gz
RUN tar -xzvf geant4-v11.2.0.tar.gz
RUN rm geant4-v11.2.0.tar.gz
WORKDIR /opt/geant4-v11.2.0
RUN mkdir build
WORKDIR /opt/geant4-v11.2.0/build
RUN apt install -y qtbase5-dev qtchooser qt5-qmake qtbase5-dev-tools
RUN cmake -DGEANT4_INSTALL_DATA=ON -DGEANT4_USE_QT=ON -DGEANT4_USE_OPENGL_X11=ON ../
RUN make -j4
RUN make install

COPY . /app
WORKDIR /app

# Make the start script executable
RUN chmod +x start_services.sh

# Expose the ports that Streamlit and API will run on
EXPOSE 8501 5000

# Start both Streamlit and the API using the script
CMD ["./start_services.sh"]