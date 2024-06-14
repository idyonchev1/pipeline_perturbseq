import os
import pandas as pd

units_table = pd.read_table("samples.txt")
samples=list(units_table.Samples.unique())

include:"rules/cellranger.smk"
include:"rules/aggr.smk"
include:"rules/multiqc.smk"

rule all:
    input:[expand("{sample}/outs/web_summary.html",sample=samples),"aggr/outs/web_summary.html","data/multiqc/multiqc_report.html"]
