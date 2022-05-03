FROM python:3.9.1-slim-buster

RUN apt-get update \ 
    && apt-get install -y cmake autoconf libtool gfortran \
                           python3-pip git gdb procps libgdal-dev \
                           libsqlite3-mod-spatialite vim valgrind curl \
    && apt-get clean -y

ENV CPLUS_INCLUDE_PATH=/usr/include/gdal
ENV C_INCLUDE_PATH=/usr/include/gdal

COPY ./requirements.txt /threedigrid-builder/

RUN pip3 install -r /threedigrid-builder/requirements.txt \
    && pip3 install GDAL==2.1.0 --global-option=build_ext \
    && pip3 install scikit-build \
    && pip3 install -r /threedigrid-builder/requirements.txt

COPY . /threedigrid-builder
RUN cd /threedigrid-builder && ./build.sh RELEASE
# Doesnt work when done from Dockerfile. Run manually!
# RUN cd /threedigrid-builder && python3 setup.py build_ext --inplace

ENV LD_LIBRARY_PATH=/usr/lib:/usr/local/lib
WORKDIR /threedigrid-builder

EXPOSE 8080