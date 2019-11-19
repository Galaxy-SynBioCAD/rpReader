FROM brsynth/rpbase

RUN apt-get install --quiet --yes --no-install-recommends \
	libxext6  \
    	libxrender-dev  && \
    conda install -y -c rdkit rdkit && \
    conda install -c conda-forge flask-restful && \

COPY rpReader.py /home/
COPY rpCache.py /home/
COPY rp2ReaderServe.py /home/

RUN python /home/rpCache.py

ENTRYPOINT ["python"]
CMD ["/home/rp2ReaderServe.py"]

# Open server port
EXPOSE 8997
