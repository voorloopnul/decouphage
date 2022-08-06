#!/usr/bin/env python3

import logging
import os
import sys
import tarfile
from pathlib import Path
import requests
from alive_progress import alive_bar

logging.basicConfig(level=20)
logger = logging.getLogger(__name__)


def download(file_url: str, file_dst: Path):
    try:

        with file_dst.open('wb') as fh_out, requests.get(file_url, stream=True) as response:
            file_length = response.headers.get('content-length')

            if file_length:
                file_length = int(file_length)
                file_length = int(file_length / 1024)

            with alive_bar(total=file_length) as bar:
                for data in response.iter_content(chunk_size=1024):
                    fh_out.write(data)
                    bar()
    except IOError:
        sys.exit(f'ERROR: Could not download file! url={file_url}, path={file_dst}')


def get_default_config_path() -> Path:
    home_path = os.path.expanduser('~')
    config_path = home_path / Path(".decouphage")
    os.makedirs(config_path, exist_ok=True)
    return config_path


def get_default_database_path(config_path) -> Path:
    database_path = config_path / Path("db")
    os.makedirs(database_path, exist_ok=True)
    return database_path


def run_download():

    default_config_path = get_default_config_path()
    default_database_path = get_default_database_path(default_config_path)
    download_file_path = default_database_path / Path("ncbi_phages.tar.gz")
    logger.info("[1/2] Downloading database...")
    download("https://labs.voorloop.com/archive/decouphage/ncbi_phages_latest.tar.gz", download_file_path)

    logger.info("[2/2] Extracting database...")
    with tarfile.open(download_file_path) as tarball:
        tarball.extractall(default_database_path)

    os.remove(download_file_path)

