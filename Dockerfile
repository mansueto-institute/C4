FROM ubuntu:16.04

RUN apt-get update
RUN apt-get install -y git ca-certificates s3cmd 
RUN apt-get_install -y libboost-all-dev libboost-doc libarmadillo-dev libarmadillo6 
RUN apt-get install -y python3 python3-pip python3-dev build-essential 
RUN apt-get install -y libgeos-dev libgdal-dev 
RUN apt-get install -y python3-gdal gdal-bin

RUN pip3 install cython matplotlib fiona pysal geopandas psycopg2
RUN git clone https://github.com/JamesSaxon/cluscious.git && cd cluscious && python3 setup.py build_ext --inplace
