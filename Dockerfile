FROM brsynth/rpcache:dev

COPY rpToolCache.py /home/
RUN python rpToolCache.py

COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY tool_rp2Reader.py /home/
COPY tool_tsvReader.py /home/
COPY tool_strReader.py /home/
