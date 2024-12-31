FROM ubuntu

RUN apt-get update && apt-get install -y \
python3 python3-pip python3-dev python3-setuptools \
g++ make libboost-all-dev wget openjdk-11-jre curl \
r-base r-base-dev


RUN curl -s https://get.nextflow.io | bash && \
mv nextflow /usr/local/bin/

# install PascalX
ADD PascalX /PascalX/
RUN mkdir -p /PascalX/build/lib
ENV DEBIAN_FRONTEND="noninteractive"
RUN echo "/PascalX/build/lib" > /etc/ld.so.conf.d/pascalx.conf
# Build
RUN cd /PascalX && make all && ldconfig && make test
RUN cd /PascalX/python/ && python3 setup.py install
# Install jupyter
RUN pip3 install jupyter

COPY . /app
WORKDIR /app

#RUN R -e "install.packages(c('WebGestaltR','optparse','stringr'),dependencies=TRUE, repos='http://cran.rstudio.com/')"

RUN pip3 install pandas statsmodels --break-system-packages
#RUN pip3 install -r requirements.txt

ENTRYPOINT ["python3", "/app/scripts/main.py", "--test"]
#ENTRYPOINT ["python3", "/app/scripts/main.py", "--trait", "maleWC", "--trait_file_path", "/app/data/maleWC.csv", "--module_directory_path", "/app/data/ker_based/"]
#ENTRYPOINT ["nextflow", "info"]
