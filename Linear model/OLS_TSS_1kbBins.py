from statsmodels.stats.multitest import multipletests
import statsmodels.formula.api as smf
import pandas as pd
import numpy as np
from tqdm import tqdm

#jobid=OLS-TSS
#sbatch --export=All -t 48:00:00 --mem=20G  -e run_log/${jobid}.err -o run_log/${jobid}.out -J $jobid --wrap="python OLS_TSS_1kbBins.py"

def accessData3(list_files, path, list_labels):
    """opens the read/transcript/gene counts of each file"""
    #create empty dictionary to store the dataframes
    dic_dataframes = {}
    #iterate to get the file names
    #make a counter to access the correct label
    count = 0
    for file in list_files:
        #open each file
        data_frame = pd.read_csv(path+file, sep='\t', index_col=0)
        #add the dataframe as an element of a dictionary with key the filename
        dic_dataframes[list_labels[count]] = data_frame
        count = count + 1

    return dic_dataframes


#access de data
#Create a variable list with the file names and one with the labels to use as dictionary keys
path='/exports/ana-geijsenlab/cmartinezmir/projects/tomoseq_timemachine/20221209_cm59_65_m002-3-4/processed/counttables/binned/'
file_names = ['M1M2M3_1000bp0ovl_cpg_coverage_table_filtered_normal_TSS3kb.bed']
labels = ['shared_filtered_merged']


#use accessData() function to obtain a dictonary with each dataset with labels as key
data_filtered = accessData3(file_names, path, labels)

#rename the header
data_filtered['shared_filtered_merged'].columns = ['start', 'end', 'Lpnp1', 'M1_s10', 'M1_s11', 'M1_s12', 'M1_s13',
       'M1_s24', 'M1_s25', 'M1_s27', 'M1_s35', 'M1_s36',
       'M1_s38', 'M2_s8', 'M2_s9', 'M2_s10', 'M2_s11', 'M2_s22', 'M2_s23',
       'M2_s24', 'M2_s25', 'M2_s36', 'M2_s37', 'M2_s38', 'M2_s39', 'M3_s1',
       'M3_s2', 'M3_s3', 'M3_s4', 'M3_s17', 'M3_s18', 'M3_s19', 'M3_s20',
       'M3_s33', 'M3_s34', 'M3_s35', 'M3_s36', 'chr_g', 'start_g', 'end_g', 'strand_g', 'gene']

data_filtered['shared_filtered_merged']=data_filtered['shared_filtered_merged'].drop(['start_g', 'end_g', 'Lpnp1', 'strand_g'], axis=1)
# merge columns 'col1', 'col2', 'col3' into a new column 'merged'
merged_data = data_filtered['shared_filtered_merged'][['chr_g', 'start', 'end', 'gene']].astype(str).agg('_'.join, axis=1)
data_filtered['shared_filtered_merged']=data_filtered['shared_filtered_merged'].drop(['start', 'end', 'chr_g', 'gene'], axis=1)
data_filtered['shared_filtered_merged']['CpG']=merged_data

# keep rows where at least one column is non-zero
df = data_filtered['shared_filtered_merged'].drop(['CpG'], axis=1)
# create mask: count how many columns are >= 4 per row, keep rows with at least 4 sections with coverage
mask = (df >= 4).sum(axis=1) >= 4
df_filtered = data_filtered['shared_filtered_merged'].loc[mask]

#build metadata
sample_metadata=pd.DataFrame(index=df_filtered.columns[:-1])
sample_metadata['Condition']=['female']*22 + ['male']*12
sample_metadata['annotation']=['central']*4+['exclude']*3+['proximal-distal']*3+ ['central']*8+['proximal-distal']*8+ ['central']*4+ ['proximal-distal']*4

# Assuming your expression table is called df_counts
# Keep only the columns with actual sample counts
sample_cols = sample_metadata.index.tolist()  # match column names to metadata

df_counts=df_filtered

# melt the table
df_long = df_counts.reset_index().melt(
    id_vars=['CpG'],  # or the column in your table that identifies the gene
    value_vars=sample_cols,
    var_name='Sample',
    value_name='Expression'
)
# sample_metadata should have 'Condition' and 'annotation'
df_long = df_long.merge(
    sample_metadata[['Condition', 'annotation']],
    left_on='Sample',
    right_index=True
)

#linear regression model
results = []
CpG = df_long['CpG'].unique()

for CpG in tqdm(CpG):
    obs_df = df_long[df_long['CpG'] == CpG].copy()

    # skip constant genes
    if np.allclose(obs_df['Expression'].values, obs_df['Expression'].values[0]):
        results.append({'CpG': CpG, 'pval': np.nan, 'note': 'constant_expression'})
        continue

    # make sure categorical
    obs_df['Condition'] = obs_df['Condition'].astype('category')
    obs_df['annotation'] = obs_df['annotation'].astype('category')
    obs_df['Sample'] = obs_df['Sample'].astype('category')

    try:
        model = smf.mixedlm(
            'Expression ~ Condition * annotation',
            data=obs_df,
            groups=obs_df['Sample'],   # group by Sample or another random effect
            re_formula="~1"
        ).fit(reml=False)

        # model.pvalues is a pandas Series
        pvals = model.pvalues

        # get Condition main effect(s)
        # keys look like 'Condition[T.male]' if your baseline is 'female'
        pval_cond = np.nan
        for term, pv in pvals.items():
            if 'Condition' in term and ':' not in term:  # exclude interaction terms
                pval_cond = pv
                break

        # get annotation main effect(s) (exclude interactions)
        pval_ann = np.nan
        for term, pv in pvals.items():
            if 'annotation' in term and ':' not in term:  # exclude interaction terms
                pval_ann = pv
                break


        # interaction term p-value
        pval_int = np.nan
        for term, pv in model.pvalues.items():
            if ':' in term and ('Condition' in term and 'annotation' in term):
                pval_int = pv
                break

    except Exception as e:
        pval_cond = np.nan
        pval_ann = np.nan
        pval_int = np.nan

    results.append({'CpG': CpG, 'pval_condition': pval_cond, 'pval_annotation': pval_ann, 'pval_interaction': pval_int})


results_df = pd.DataFrame(results)
# Save to Excel
results_df.to_excel(path+"OLS_TSS.xlsx")
