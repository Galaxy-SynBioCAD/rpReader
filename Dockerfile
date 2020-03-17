FROM brsynth/rpcache

COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY rpToolCache.py /home/
COPY tool_rp2Reader.py /home/

RUN python rpToolCache.py
