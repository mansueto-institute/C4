FROM ubuntu:16.04

RUN apt-get update
RUN apt-get install -y zip git ca-certificates s3cmd
RUN apt-get install -y libboost-all-dev libboost-doc libarmadillo-dev libarmadillo6
RUN apt-get install -y python3 python3-pip python3-dev build-essential
RUN apt-get install -y libgeos-dev libgdal-dev
RUN apt-get install -y python3-gdal gdal-bin

RUN pip3 install cython matplotlib fiona pysal geopandas psycopg2

# https://github.com/settings/tokens/
RUN mkdir C4 && cd C4 && git init && git pull https://MY_PERSONAL_TOKEN@github.com/JamesSaxon/C4.git && python3 setup.py build_ext --inplace

RUN touch ~/.netrc

ENV AWS_DEFAULT_REGION=us-east-1 \
    AWS_ACCESS_KEY_ID=MY_AWS_ACCESS_KEY_ID \
    AWS_SECRET_ACCESS_KEY=MY_AWS_SECRET_ACCESS_KEY

CMD cd C4 && \
    echo RUNNING :: $STATE $SEED $METHOD && \
    ./run_iter.sh 2>&1 | tee ${STATE}-${SEED}-${METHOD}.out && \
    zip -r $(printf "%s_s%03d%s.zip" $STATE $SEED $([ "$METHOD" != "" ] && echo _$METHOD)) \
           ${STATE}-${SEED}-${METHOD}.out res/ && \
    s3cmd put *zip s3://jsaxon-test-bucket/
