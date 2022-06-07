FROM ubuntu:20.04
ENV TZ=Europe/Copenhagen
ENV DEBIAN_FRONTEND=noninteractive

RUN apt update && apt -y install wget unzip nano tree python3 python3-pip python3-distutils python-is-python3 \
    r-base tree prodigal ncbi-blast+ && \
    rm -rf /var/lib/apt/lists/* && apt clean

RUN mkdir /data

COPY requirements.txt /data/
RUN pip install -r /data/requirements.txt

COPY decouphage.py /data/
COPY decouphage_db.py /data/
COPY src /data/src

WORKDIR data
RUN python decouphage_db.py download
ENTRYPOINT ["python3", "decouphage.py"]
CMD ["--help"]
