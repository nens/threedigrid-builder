FROM python:3.6-slim-stretch as build

RUN apt-get update \ 
    && apt-get install -y cmake autoconf libtool gfortran libnetcdf-dev libnetcdff-dev libgdal-dev libproj-dev libshp-dev python3-pip git\
    && apt-get clean -y

RUN pip3 install numpy ipython jupyter

RUN git clone --branch v2.6 --depth 1 https://github.com/dcesari/fortrangis.git \
    && cd /fortrangis && autoreconf -vi \
    && ./configure --disable-doxydoc \
    && make install

COPY ./threedi-calculationcore /threedicore
RUN cd /threedicore && ./full_build.sh

COPY ./py3di /py3di
RUN cd /py3di && python3 setup.py install



FROM python:3.6-slim-stretch as deps

RUN apt-get update && apt-get install libnetcdf11 libnetcdff6 libgdal20 libshp2 liblas3 liblapack3 libgomp1 -y



FROM python:3.6-slim-stretch

COPY --from=build /usr/local/lib/libfortran*.so* /usr/local/lib/ 
COPY --from=build /opt/threedicore/ /opt/threedicore/
COPY --from=build /usr/local/lib/python3.6/site-packages/ /usr/local/lib/python3.6/site-packages/
COPY --from=build /usr/local/bin/ipython /usr/local/bin/ipython
COPY --from=build /usr/local/bin/*jupyter* /usr/local/bin/
COPY --from=build /usr/local/etc/*jupyter* /usr/local/etc/
COPY --from=build /usr/local/share/*jupyter* /usr/local/share/
COPY --from=deps /lib/x86_64-linux-gnu/*.so* /lib/x86_64-linux-gnu/
COPY --from=deps /etc/alternatives/libblas* /etc/alternatives/
COPY --from=deps /etc/alternatives/liblapack* /etc/alternatives/
COPY --from=deps /usr/lib/x86_64-linux-gnu/*.so* /usr/lib/x86_64-linux-gnu/
COPY --from=deps /usr/lib/libarmadillo*.so* /usr/lib/
COPY --from=deps /usr/lib/libogdi*.so* /usr/lib/
COPY --from=deps /usr/lib/libgdal*.so* /usr/lib/
COPY --from=deps /usr/share/gdal/2.1/* /usr/share/gdal/2.1/
COPY --from=deps /usr/lib/libblas* /usr/lib/
COPY --from=deps /usr/lib/liblapack* /usr/lib/
COPY --from=deps /usr/lib/libmfhdfalt*.so* /usr/lib/
COPY --from=deps /usr/lib/libdfalt*.so* /usr/lib/
COPY --from=deps /usr/lib/libarpack*.so* /usr/lib/


ENV LD_LIBRARY_PATH=/usr/lib:/usr/local/lib:/opt/threedicore/lib
#ENTRYPOINT ["jupyter", "notebook", "--ip=0.0.0.0", "--allow-root"]
