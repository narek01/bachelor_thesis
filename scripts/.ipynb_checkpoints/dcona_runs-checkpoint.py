import multiprocessing
import gspread_pandas
import telegram_send
import pandas as pd
import numpy as np
import pathlib
import dcona
from datetime import datetime
import time

class Experiment():
    def __init__(self, *args, **kwargs):
        self._init_args = args
        self._init_kwargs = kwargs
        self.gspread = MetaExporter()
        self._run_number = (int(self.gspread._df_cache.index.max()) + 1 if len(self.gspread._df_cache) else 1)
        self._meta_row = self.gspread.create_row(kwargs)
        self.replace_in_row('status', 'New')
        
    def ztest_wrapped(self, *args, **kwargs):
        self._ztest_started = datetime.now()
        self.replace_in_row('ztest_start', self._ztest_started)
        telegram_send.send(messages=[f'Ztest started, experiment {self._init_kwargs["experiment"]}'])
        try:
            args_to_export = ['reference_group', 'experimental_group']
            self.replace_in_row('ztest_args', {key: kwargs['ztest_kwargs'][key] for key in args_to_export})
            self.replace_in_row('status', 'Running')
            self.ztest_result = dcona.ztest(**kwargs['ztest_kwargs'])
            self._ztest_finished = datetime.now()
            self._ztest_time_spent = (self._ztest_finished) - (self._ztest_started)
            dirpath = pathlib.Path(f'/home/jovyan/shared/narek/outputs/experiments/{self._init_kwargs["experiment"]}')
            dirpath.mkdir(parents=True, exist_ok=True)
            filepath = dirpath/kwargs["filename"]
            self.ztest_result.to_csv(filepath)
            telegram_send.send(messages=[f'Ztest performed (exp {self._init_kwargs["experiment"]}), time spent {self._ztest_time_spent}'])
            self.replace_in_row('ztest_path', filepath)
            self.replace_in_row('ztest_end', self._ztest_finished)
            self.replace_in_row('ztest_time', self._ztest_time_spent)
            self.replace_in_row('status', 'Finished')
        except Exception as exc:
            self.replace_in_row('status', 'Error')
            telegram_send.send(messages=[f'Error occured:\n {exc}'])
            raise Exception
        
    def run_ztest(self, *args, **kwargs):
        multiprocessing.Process(target=self.ztest_wrapped, args=args, kwargs=kwargs).start()
        
    def replace_in_row(self, column, value):
        self._meta_row.iloc[0].loc[column] = value
        self.gspread.send_row(self._meta_row, self._run_number)

class DataProcessor:
    def __init__(self, expression_df: pd.DataFrame, annotation_df: pd.DataFrame, mirna_df: pd.DataFrame):
        self.annotation = annotation_df.copy()
        self.gene_exp = expression_df.copy()
        self.mirna = mirna_df.copy()
    
    def cutoff_expressions(self, subtypes: list = None, cutoff: float = 95, remove_tr_version: bool = True):
        if subtypes:
            ids = self.annotation.loc[self.annotation['Sample type'].isin(subtypes)].index.to_numpy()
            exp_all = extract_samples(self.gene_exp, ids, additional_columns=['gene symbol'])
        else:
            exp_all = self.gene_exp.copy()
        self.cutted_exp = filter_by_median(exp_all, additional_columns=['gene symbol'], cutoff=cutoff)
        if remove_tr_version:
            self.cutted_exp = remove_transcript_version(self.cutted_exp)
        print(self.cutted_exp.shape[0], 'genes and', self.cutted_exp.shape[1], 'samples are chosen')

    def cutoff_mirnas(self, cutoff: float = 99):
        mirna_5prime = isomir_groupby_5prime(self.mirna)
        self.mirna_5prime = filter_by_median(mirna_5prime, additional_columns=['median'], cutoff=cutoff)
        print(self.mirna_5prime.shape[0], 'miRNAs are chosen')
        
    def final_data(self, description_df: pd.DataFrame, raw_interaction_df: pd.DataFrame):
        self.expr_united = pd.concat(tcga_match_samples(self.cutted_exp.drop('gene symbol', axis=1), self.mirna_5prime))
        self.interaction_df = raw_interaction_df[['isomir', 'transcript']].drop_duplicates()
        self.interaction_df.columns = ['Source', 'Target']
        return brush_data_before_dcona(self.expr_united, description_df, self.interaction_df)
        
class CuttingSamples:
    def __init__(self, exp_df, extra_columns=None):
        self._exp_df = exp_df
        self._extra_columns = extra_columns
    
    def cut(self, transcript, q=10, plot=True, suppress_print=False):
        self.__transcript = transcript
        df = self._exp_df.loc[[transcript]].drop(self._extra_columns, axis=1)
        series = df.T.squeeze()
        self._qcut_df = pd.qcut(series, q=q, duplicates='drop')
        if plot:
            df.T.plot.hist(bins=20, alpha=0.5)
        if not suppress_print:
            print(self._qcut_df.value_counts(sort=False))
        
    def extract_samples(self, left_num=1, right_num=1, return_values=False, return_all=False):
        left_categories = self._qcut_df.cat.categories[:left_num]
        right_categories = self._qcut_df.cat.categories[-right_num:]
        left_samples = self._qcut_df.loc[self._qcut_df.isin(left_categories)]
        right_samples = self._qcut_df.loc[self._qcut_df.isin(right_categories)]
        left_df = self._df_sample_group(left_samples, "Low_expr")
        right_df = self._df_sample_group(right_samples, "High_expr")
        self.__cutted_samples = pd.concat([left_df, right_df]).reset_index(drop=True)
        if return_all:
            other_samples = self._exp_df.loc[:, ~self._exp_df.columns.isin(self.__cutted_samples.Sample)].drop(self._extra_columns, axis=1).T
            other_df = self._df_sample_group(other_samples, "_Other")
            self.__cutted_samples = pd.concat([self.__cutted_samples, other_df]).reset_index(drop=True)
        if return_values:
            self.__cutted_samples['Value'] = self._exp_df.loc[self.__transcript, self.__cutted_samples['Sample']].to_numpy()
        return self.__cutted_samples
    

    def _df_sample_group(self, samples, group):
        df = pd.DataFrame({
            "Sample": list(samples.index),
            "Group": [group]*len(samples)
        })
        return df

    def save_annotation(self, gene_name, directory="/home/jovyan/shared/narek/outputs/sample_cuts"):
        filename = f"{gene_name}_{self.__transcript}_{len(self.__cutted_samples)}_samples.csv"
        path = pathlib.Path(directory)/pathlib.Path(filename)
        if path.exists():
            print('File already exists:', path)
            return
        self.__cutted_samples.to_csv(path)
        print('Annotation is saved at:', path)


class MetaExporter:
    def __init__(self, account_file='/home/jovyan/diploma_scripts/.service_account.json'):
        my_config = gspread_pandas.conf.get_config('/'.join(account_file.split('/')[:-1]), account_file.split('/')[-1])
        self.spread = gspread_pandas.Spread('DCoNA_runs', config=my_config)
        self.df
    
    @property
    def df(self):
        df = self.spread.sheet_to_df()
        df.index = df.index.astype(int)
        self._df_cache = df
        return df
    
    def create_row(self, dct={}):
        row = pd.DataFrame(data=dct, index=[0], columns=self._df_cache.columns)
        return row
    
    def send_row(self, row, row_number):
        df = self.df
        row.index = [row_number]
        if row_number in df.index:
            df.loc[row_number] = row.to_numpy()[0]
        else:
            df = pd.concat([df, row], ignore_index=False)
        df = df.sort_index()
        self.spread.df_to_sheet(df, freeze_headers=True)
        return df
      

### Utils
def extract_samples(df, sample_ids, additional_columns=None):
    if additional_columns:
        return df.loc[:, np.concatenate((additional_columns, sample_ids))]
    return df.loc[:, sample_ids]

def filter_by_median(df, additional_columns=None, max_rows=None, cutoff=None):
    filtering_df = df.copy(); df_out = df.copy() 
    if additional_columns:
        filtering_df = filtering_df.drop(additional_columns, axis=1)
    df_out['median'] = filtering_df.median(axis=1)
    df_out = df_out.loc[df_out['median'] > 0]
    df_out.sort_values(['median'], ascending=False, inplace=True)
    if cutoff:
        temp_df = df_out.drop(additional_columns, axis=1)
        temp_df = 2**temp_df - 1
        cutoff_value = temp_df.median(axis=1).sort_values(ascending=False).cumsum()[-1] * float(cutoff) / 100
        df_out = df_out.loc[temp_df.median(axis=1).sort_values(ascending=False).cumsum() < cutoff_value]
    if max_rows:
        df_out = df_out.iloc[:max_rows, :]
    try:
        df_out.drop('median', axis=1, inplace=True)
    except: pass
    return df_out

def isomir_groupby_5prime(isomir, counts_data=False):
    isomir["median"] = isomir.iloc[:, 3:].median(axis=1)
    isomir = isomir.sort_values("median", ascending=False)
    # isomir = isomir.loc[isomir["median"] > 0]
    isomirs_with_mirdb = set(open("/home/jovyan/shared/miRNA_predictions_BRCA/isomiRs_with_miRDB.txt").read().strip().split())
    if not counts_data:
        isomir.iloc[:, 3:] = 2**isomir.iloc[:, 3:] - 1
    mature_5p = []
    for mm in isomir["mature"]:
        if type(mm) != str:
            mature_5p.append(None)
            continue

        has_mirdb = False
        for m in mm.split(", "):
            m_5p = "|".join(m.split("|")[:-1])
            if m_5p in isomirs_with_mirdb:
                has_mirdb = True
                mature_5p.append(m_5p)
                break

        if not has_mirdb:
            mature_5p.append(None)

    isomir["mature"] = mature_5p
    isomir = isomir.loc[isomir["mature"].notna()]
    isomir = isomir.drop(columns=["hairpin", "repeat", "median"])
    isomir = isomir.groupby("mature").sum()
    if not counts_data:
        isomir = np.log2(isomir + 1)
    isomir["median"] = isomir.median(axis=1)
    isomir = isomir.loc[isomir["median"] > 0]
    isomir = isomir.sort_values("median", ascending=False)

    return isomir

def brush_data_before_dcona(expr, annotation, interaction):
    expr, annotation, interaction = expr.copy(), annotation.copy(), interaction.copy()
    common_samples = list(set(expr.columns) & set(annotation["Sample"]))
    common_genes = list((
            set(interaction['Source'].unique()) | 
            set(interaction['Target'].unique())
        ) & set(expr.index))
    common_mirnas = [c for c in common_genes if c.startswith("hsa")]
    common_transcripts = [c for c in common_genes if c.startswith("ENS")]
    interaction = interaction.loc[interaction['Source'].isin(common_mirnas) & interaction['Target'].isin(common_transcripts)]
    annotation = annotation.loc[annotation['Sample'].isin(common_samples)]
    expr = expr.loc[common_genes, common_samples]
    return expr, annotation, interaction

def remove_transcript_version(expr):
    expr.index = expr.index.map(lambda x: x.split('.')[0])
    return expr

def tcga_match_samples(df1, df2):
    tech_cols1 = [c for c in df1.columns if not c.startswith("TCGA")]
    tech_cols2 = [c for c in df2.columns if not c.startswith("TCGA")]
    common_samples = sorted(list(set(df1.columns) & set(df2.columns)))
    return df1[tech_cols1 + common_samples], df2[tech_cols2 + common_samples]

def sample_df_for_swarplot(expr_dict, transcript, gene_name, subtype_dict, extra_columns=['gene symbol']):
    output_df_list = list()
    for subtype in subtype_dict.keys():
        cutter = CuttingSamples(expr_dict[subtype], extra_columns=extra_columns)
        cutter.cut(transcript, q=subtype_dict[subtype]['qcut'], plot=False, suppress_print=True)
        df = cutter.extract_samples(left_num=subtype_dict[subtype]['left'], 
                                    right_num=subtype_dict[subtype]['right'], 
                                    return_values=True, return_all=True)
        df['Subtype'] = [subtype] * df.shape[0]
        output_df_list.append(df)
    return pd.concat(output_df_list)