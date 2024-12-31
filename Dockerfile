FROM ubuntu:20.04

RUN apt-get update
# install python dependencies
RUN apt-get install -y --no-install-recommends \
python3 python3-pip python3-dev python3-setuptools openjdk-11-jre curl

# install R dependencies
#RUN apt-get install -y \
#r-base r-base-dev


# build PascalX
COPY ./PascalX /PascalX
RUN mkdir -p /PascalX/build/lib
RUN export DEBIAN_FRONTEND=noninteractive
ARG DEBIAN_FRONTEND=noninteractive
ENV TZ=America/Chicago
# install PascalX dependencies
RUN apt-get update && apt-get install -y --no-install-recommends \
g++ make libboost-all-dev wget unzip \
libjpeg-dev zlib1g-dev libfreetype6-dev liblcms2-dev \
libopenjp2-7-dev libtiff5-dev libwebp-dev libimagequant-dev libxcb1-dev
RUN echo "/PascalX/build/lib" > /etc/ld.so.conf.d/pascalx.conf
RUN pip3 install "pillow>=6.2.0,<11.0" "matplotlib==3.5.0" "pybind11" "numpy==1.19.2"
# Build
RUN cd /PascalX && make all && ldconfig && make test
RUN cd /PascalX/python/ && python3 setup.py install
RUN pip3 install jupyter

CMD ["python3","-c","'from PascalX import genescorer'"]


# install nextflow
RUN curl -s https://get.nextflow.io | bash && \
mv nextflow /usr/local/bin/



ENTRYPOINT ["python3", "/app/scripts/main.py", "--test"]
#ENTRYPOINT ["python3", "/app/scripts/main.py", "--trait", "maleWC", "--trait_file_path", "/app/data/maleWC.csv", "--module_directory_path", "/app/data/ker_based/"]
#ENTRYPOINT ["nextflow", "info"]
