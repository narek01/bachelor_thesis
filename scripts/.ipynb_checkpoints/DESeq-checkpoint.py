import pandas as pd
import rpy2.robjects as ro
from rpy2.robjects.packages import importr
from rpy2.robjects import pandas2ri
from rpy2.robjects.conversion import localconverter

pandas2ri.activate()
DESeq2 = importr("DESeq2")

from rpy2.robjects import default_converter
from rpy2.robjects.conversion import rpy2py
base = importr("base")

from pathlib import Path
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from time import time, gmtime, sleep
import sys
sys.path.append("/home/jovyan/HSE-Bioinformatics")
sys.path.append("/home/jovyan/diploma_scripts/scripts")
from bio import *
import dcona_runs
import dcona
from tqdm.auto import tqdm; tqdm.pandas()
import telegram_send


def deseq(meta: pd.DataFrame, counts: pd.DataFrame, formula: str, ref: str, exp: str):
    dds = DESeq2.DESeqDataSetFromMatrix(countData=counts, colData=meta, design=ro.Formula(formula))
    dds = DESeq2.DESeq(dds)
    print(f"Group_{exp}_vs_{ref}")
    resR = DESeq2.results(dds, name=f"Group_{ref}_vs_{exp}")
    res = r_to_df(resR)
    res = res.sort_values("padj")
    res = res.loc[res["padj"] < 0.05]
    res = res.loc[res["log2FoldChange"].abs() > 0.5]

    return res


def r_to_df(r_df):
    with localconverter(default_converter + pandas2ri.converter):
        return rpy2py(base.as_data_frame(r_df))


SHARED_PATH = Path('/home/jovyan/shared')
TCGA_PATH = SHARED_PATH / 'TCGA_data'
OUTPUT_PATH = SHARED_PATH / 'narek/outputs'
BRCA_DATA = SHARED_PATH / 'narek/Diplom2/data_BRCA'

TCGA_STUDY = 'TCGA-BRCA'
iso_type = "pan_cancer_exclusive_log2_FPM_DESeq2"
annotation = rt(BRCA_DATA/'annotation.tsv')
gene_ids = rc('~/diploma_scripts/aux/expr95_mirna99.csv').index

gene_exp_counts = rt(f"{TCGA_PATH}/{TCGA_STUDY}/RSEM_transcript_pan_cancer_counts.tsv")
mirna_counts = rt(f"{TCGA_PATH}/{TCGA_STUDY}/isoMiRmap_pan_cancer_exclusive_counts.tsv")
mirna_counts_5prime = dcona_runs.isomir_groupby_5prime(mirna_counts, counts_data=True)
united_data = pd.concat(
    tcga_match_samples(
        dcona_runs.remove_transcript_version(gene_exp_counts.drop('gene symbol', axis=1)),
        mirna_counts_5prime.drop('median', axis=1)
    )).loc[gene_ids]

del gene_exp_counts; del mirna_counts; del mirna_counts_5prime

distrib_dict = {
    '1': 'AGO2_normal_ENST00000220592_62_samples.csv',
    '2': 'AGO2_luminalA_ENST00000220592_113_samples.csv',
    '3': 'AGO2_luminalB_ENST00000220592_84_samples.csv',
    '4': 'AGO2_basalLike_ENST00000220592_74_samples.csv',
    '5': 'DROSHA_normal_ENST00000513349_62_samples.csv',
    '6': 'DROSHA_luminalA_ENST00000513349_150_samples.csv',
    '7': 'DICER1_normal_ENST00000393063_62_samples.csv',
    '8': 'DICER1_luminalA_ENST00000393063_150_samples.csv',
    '9': 'DGCR8_luminalA_ENST00000351989_150_samples.csv',
    '10': 'DGCR8_normal_ENST00000351989_62_samples.csv'
}

for key in distrib_dict.keys():
    descr = rc(SHARED_PATH/f'narek/outputs/sample_cuts/{distrib_dict.get(key)}')
    descr = descr.set_index('Sample')
    extracted_united = dcona_runs.extract_samples(united_data, descr.index)
    res = deseq(meta=descr, counts=extracted_united, formula="~Group", ref="Low_expr", exp="High_expr")
    res.to_csv(OUTPUT_PATH/'experiments'/key/'DESeq2.tsv', sep='\t')