import os
import pandas as pd
from pathlib import Path

paths = [
    "qc",
    "qc/{sample}",
    "rawdata/{sample}",
    "logs/{sample}",
    "outs/{sample}/trimmedReads",
    "outs/{sample}/{counttype}",
    "outs/agg/{counttype}",
    "results/{sample}/plots",
    "results/agg/clusterprofiler/{contrast}",
    "results/agg/deseq2/{counttype}/{contrast}",
    "results/agg/deseq_telescope/{tecounttype}/{contrast}",
    "results/agg/repeatanalysis_telescope"
    ]
paths = paths + []
for path in paths:
    for sample in samples:
        for counttype in config["counttypes"]:
            for tecounttype in config["tecounttypes"]:
                for contrast in config["contrasts"]:
                    os.makedirs(path.format(sample=sample, counttype=counttype,tecounttype=tecounttype, contrast=contrast), exist_ok = True)
