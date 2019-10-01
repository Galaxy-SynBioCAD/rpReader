FROM ibisba/rpbase

RUN apt-get install --quiet --yes --no-install-recommends \ 
	libxext6  \
    	libxrender-dev  && \
    conda install -y -c rdkit rdkit && \
    mkdir input_cache && \
    wget https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_xref.tsv -P /home/input_cache/ && \
    wget https://www.metanetx.org/cgi-bin/mnxget/mnxref/reac_xref.tsv -P /home/input_cache/ && \
    wget https://www.metanetx.org/cgi-bin/mnxget/mnxref/chem_prop.tsv -P /home/input_cache/ && \
    wget https://www.metanetx.org/cgi-bin/mnxget/mnxref/comp_xref.tsv -P /home/input_cache/

COPY rpReader.py /home/

#get the rules
RUN wget https://retrorules.org/dl/preparsed/rr02/rp3/hs -O /home/rules_rall_rp3.tar.gz && \
    tar xf /home/rules_rall_rp3.tar.gz -C /home/ && \
    mv /home/retrorules_rr02_rp3_hs/retrorules_rr02_flat_all.tsv /home/input_cache/rules_rall.tsv && \
    rm -r /home/retrorules_rr02_rp3_hs && \
    rm /home/rules_rall_rp3.tar.gz

#get the additional information
RUN wget https://retrorules.org/dl/this/is/not/a/secret/path/rr02 -O /home/rr02_more_data.tar.gz && \
    tar xf /home/rr02_more_data.tar.gz -C /home/ && \
    mv /home/rr02_more_data/compounds.tsv /home/input_cache/rr_compounds.tsv && \
    rm -r /home/rr02_more_data && \
    rm /home/rr02_more_data.tar.gz

RUN python rpReader.py
