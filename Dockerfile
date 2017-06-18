FROM ubuntu:16.04

RUN apt-get update
RUN apt-get install -y zip git ca-certificates s3cmd 
RUN apt-get install -y libboost-all-dev libboost-doc libarmadillo-dev libarmadillo6 
RUN apt-get install -y python3 python3-pip python3-dev build-essential 
RUN apt-get install -y libgeos-dev libgdal-dev 
RUN apt-get install -y python3-gdal gdal-bin

RUN pip3 install cython matplotlib fiona pysal geopandas psycopg2
RUN git clone https://github.com/JamesSaxon/cluscious.git && cd cluscious && python3 setup.py build_ext --inplace

RUN touch ~/.netrc

ENV AWS_DEFAULT_REGION=us-east-1 \
    AWS_ACCESS_KEY_ID=YOUR_KEY_HERE \
    AWS_SECRET_ACCESS_KEY=YOUR_SECRET_KEY_HERE

CMD cd cluscious && \
    ./run_iter.sh $STATE $SEED && \
    zip -r $(printf "%s_s%03d.zip" $STATE $SEED) res/ && \
    s3cmd put *zip s3://jsaxon-test-bucket/

