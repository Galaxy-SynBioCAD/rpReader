FROM brsynth/rpbase

RUN conda install -c conda-forge flask-restful

COPY rpReader.py /home/
COPY rp2ReaderServe.py /home/

ENTRYPOINT ["python"]
CMD ["/home/rp2ReaderServe.py"]

# Open server port
EXPOSE 8997
