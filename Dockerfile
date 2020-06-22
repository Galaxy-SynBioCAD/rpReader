FROM brsynth/rpcache:dev

COPY rpTool.py /home/
COPY rpToolServe.py /home/
COPY galaxy/tool_rp2Reader.py /home/
COPY galaxy/tool_tsvReader.py /home/
#COPY tool_strReader.py /home/
