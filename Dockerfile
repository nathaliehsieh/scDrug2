FROM python:3.8-slim-buster
WORKDIR single-cell-analysis
COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt
COPY . .
RUN apt-get update && apt-get install -y git unzip wget

## scMatch
RUN git clone https://github.com/asrhou/scMatch.git /opt/scMatch

RUN unzip /opt/scMatch/refDB/FANTOM5/10090_HID.csv.zip -d /opt/scMatch/refDB/FANTOM5/ && rm /opt/scMatch/refDB/FANTOM5/10090_HID.csv.zip
RUN unzip /opt/scMatch/refDB/FANTOM5/10090_symbol.csv.zip -d /opt/scMatch/refDB/FANTOM5/ && rm /opt/scMatch/refDB/FANTOM5/10090_symbol.csv.zip
RUN unzip /opt/scMatch/refDB/FANTOM5/9606_HID.csv.zip -d /opt/scMatch/refDB/FANTOM5/ && rm /opt/scMatch/refDB/FANTOM5/9606_HID.csv.zip
RUN unzip /opt/scMatch/refDB/FANTOM5/9606_symbol.csv.zip -d /opt/scMatch/refDB/FANTOM5/ && rm /opt/scMatch/refDB/FANTOM5/9606_symbol.csv.zip 

# RUN wget -P /opt/scMatch/refDB/ https://ndownloader.figshare.com/articles/7397003?private_link=efd2969ce20fae5c118f
# RUN mv /opt/scMatch/refDB/7397003?private_link=efd2969ce20fae5c118f /opt/scMatch/refDB/7397003.zip
# RUN unzip /opt/scMatch/refDB/7397003.zip -d /opt/scMatch/refDB/FANTOM5_SingleR
# RUN rm -rf /opt/scMatch/refDB/7397003.zip

# RUN cp /opt/scMatch/refDB/FANTOM5/human_samples_oto.txt /opt/scMatch/refDB/FANTOM5/mouse_samples_oto.txt /opt/scMatch/refDB/FANTOM5/hgmm_samples_oto.txt /opt/scMatch/refDB/FANTOM5_SingleR/

RUN sed -i 's/\.ix/.loc/g' /opt/scMatch/scMatch.py
RUN sed -i 's/loc\[commonRows, ].fillna(0\.0)/reindex(commonRows, axis="index", fill_value=0.0)/g' /opt/scMatch/scMatch.py

## CaDRReS-Sc
RUN git clone https://github.com/CSB5/CaDRReS-Sc.git /opt/CaDRReS-Sc
RUN wget https://www.dropbox.com/s/3v576mspw5yewbm/GDSC_exp.tsv -O /opt/CaDRReS-Sc/data/GDSC/GDSC_exp.tsv

RUN sed -i 's/import tensorflow as tf/import tensorflow.compat.v1 as tf\ntf.disable_v2_behavior()/g' /opt/CaDRReS-Sc/cadrres_sc/model.py
RUN sed -i 's/import tensorflow\.python\.util\.deprecation as deprecation/from tensorflow.python.util import deprecation/g' /opt/CaDRReS-Sc/cadrres_sc/model.py

CMD [ "/bin/bash" ]


