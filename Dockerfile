FROM python:3.8-slim-buster
WORKDIR single-cell-analysis
COPY requirements.txt requirements.txt
RUN pip3 install -r requirements.txt
COPY . .

## scMatch installation
RUN apt-get update && apt-get install -y git unzip wget
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

CMD [ "/bin/bash" ]


