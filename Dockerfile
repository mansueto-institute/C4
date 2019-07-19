FROM ubuntu:16.04

RUN apt-get update
RUN apt-get install -y zip git ca-certificates vim
RUN apt-get install -y libboost-all-dev libboost-doc libarmadillo-dev libarmadillo6
RUN apt-get install -y python3 python3-pip python3-dev build-essential
RUN apt-get install -y libgeos-dev libgdal-dev
RUN apt-get install -y python3-gdal gdal-bin

RUN pip3 install cython==0.28.5 matplotlib==2.2.2 fiona==1.8.0 pysal==1.14.4 geopandas==0.4.0 psycopg2==2.6.2 pandas==0.23.4 scipy==1.1.0 pyproj==1.9.5.1 descartes

RUN git clone https://github.com/JamesSaxon/C4.git && cd C4 && python3 setup.py build_ext --inplace

RUN touch ~/.netrc

CMD cd C4 && ./run_iter.sh 

