FROM ubuntu:22.04
ENV TZ=Europe/Copenhagen
ENV DEBIAN_FRONTEND=noninteractive

RUN apt update && apt -y install python3 python3-pip python-is-python3 prodigal ncbi-blast+ trnascan-se less && \
    rm -rf /var/lib/apt/lists/* && apt clean

RUN mkdir /data

COPY requirements.txt /data/
RUN pip install -r /data/requirements.txt

COPY decouphage /data/
COPY src /data/src

WORKDIR data
RUN python decouphage --download_db
CMD ["python", "decouphage.py"]
