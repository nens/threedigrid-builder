FROM python:3.9.1-slim-buster as build

RUN apt-get update \ 
    && apt-get install -y cmake autoconf libtool gfortran \
                          libgdal-dev python3-pip git procps \
    && apt-get clean -y

ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
ENV C_INCLUDE_PATH=/usr/include/gdal

COPY ./requirements.txt /threedicore/

RUN pip3 install -r /threedicore/requirements.txt \
    && pip3 install GDAL==2.1.0 --global-option=build_ext \
    && pip3 install cython

COPY . /gridgenerator
RUN cd /gridgenerator && ./full_build.sh RELEASE
#RUN cd /gridgenerator/py3digrd && python3 setup.py


FROM python:3.9.1-slim-buster

LABEL name=gridgenerator3di

ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
ENV C_INCLUDE_PATH=/usr/include/gdal

COPY ./requirements.txt /tmp/

RUN apt-get update && apt-get install gfortran python3-pip libgomp1 -y \
    && pip3 install -r /tmp/requirements.txt \
    && apt-get remove --purge build-essential python3-pip -y \
    && apt-get autoremove --purge -y \
    && rm -rf /var/lib/apt/lists/*

COPY --from=build /usr/local/lib/libgridgen*.so* /usr/local/lib/ 
# Copy GDAL python files
COPY --from=build /usr/local/lib/python3.7/site-packages/GDAL-2.1.0-py3.7.egg-info \
    /usr/local/lib/python3.7/site-packages/GDAL-2.1.0-py3.7.egg-info
COPY --from=build /usr/local/lib/python3.7/site-packages/gdal* /usr/local/lib/python3.7/site-packages/
COPY --from=build /usr/local/lib/python3.7/site-packages/ogr.py /usr/local/lib/python3.7/site-packages/
COPY --from=build /usr/local/lib/python3.7/site-packages/osr.py /usr/local/lib/python3.7/site-packages/
COPY --from=build /usr/local/lib/python3.7/site-packages/osgeo /usr/local/lib/python3.7/site-packages/osgeo

ENV LD_LIBRARY_PATH=/usr/lib:/usr/local/lib

