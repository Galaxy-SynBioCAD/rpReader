#FROM brsynth/rpcache-rest
FROM brsynth/rpcache

COPY rpTool.py /home/
COPY rpToolServe.py /home/
RUN python rpToolCache.py
