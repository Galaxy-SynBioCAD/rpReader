FROM brsynth/rpcache:v2

COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY galaxy/code/tool_rp2Reader.py /home/
COPY galaxy/code/tool_tsvReader.py /home/
