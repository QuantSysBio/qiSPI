shell.executable("/bin/bash")

import pandas as pd
import os, multiprocessing
import yaml
from snakemake import utils
from snakemake.utils import min_version
from snakemake import logging

min_version("6.0")
shell.prefix("set -euo pipefail;")

config = yaml.load(open("INPUT/config.yaml", "r+"), Loader=yaml.FullLoader)
dependencies = yaml.load(open("SOURCE/dependencies.yaml", "r+"), Loader=yaml.FullLoader)

snakefiles = "SOURCE/"
include: snakefiles + "rules.py"

rule all:
    input:
        split_quantitation = expand("OUTPUT/tmp/{protein_name}/split_quantitation.txt",
            protein_name=config["protein_name"]),
        finalfilteredResults = expand("OUTPUT/{protein_name}/filteredResults_final.RData",
            protein_name=config["protein_name"]),
        finalfilteredMeans = expand("OUTPUT/{protein_name}/filteredMeans_final.RData",
            protein_name=config["protein_name"]),
        KineticsDB = "OUTPUT/KineticsDB.csv"