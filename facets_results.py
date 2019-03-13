import pandas as pd
import numpy as np
import seaborn as sns
import os
import re


desired_width = 320
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width)
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.expand_frame_repr', False)


fname = 'facets_fga_all_20190311.csv'
path = '/Users/fongc2/Documents/github/MSK/facets-suite/'
pathfilename_folders = os.path.join(path, fname)
df_facets = pd.read_csv(pathfilename_folders, header=0, low_memory=False)

# Load portal fga data
fname = 'Mutation_Count_vs_Fraction_of_Genome_Altered.tsv'
path = '/Users/fongc2/Documents/github/MSK/mskimpact'
pathfilename_folders = os.path.join(path, fname)
df_cna_portal = pd.read_csv(pathfilename_folders, header=0, sep='\t')

# Get sample ids
s_id = df_facets['FOLDER_NAME_TOP'].apply(lambda x: x[:17])
df_facets = df_facets.assign(SAMPLE_ID=s_id)
# Clean by replacing Inf with NaN
df_facets = df_facets.replace([np.Inf, -np.Inf], np.NaN)

# Add column converting purity into 3 categories
df_facets = df_facets.assign(PURITY_LEVEL=df_facets['purity'])
purity_levels = [0.3, 0.8]
purity_low = df_facets['PURITY_LEVEL'] < purity_levels[0]
purity_good = (df_facets['PURITY_LEVEL'] >= purity_levels[0]) & (df_facets['PURITY_LEVEL'] < purity_levels[1])
purity_high = df_facets['PURITY_LEVEL'] > purity_levels[1]
df_facets.loc[purity_low, 'PURITY_LEVEL'] = '<%s' % purity_levels[0]
df_facets.loc[purity_good, 'PURITY_LEVEL'] = '>=%s & <%s' % (purity_levels[0], purity_levels[1])
df_facets.loc[purity_high, 'PURITY_LEVEL'] = '>%s' % purity_levels[1]
# Label unknown purity as 'Not Calculated'
df_facets.loc[df_facets['purity'].isnull(), 'PURITY_LEVEL'] = 'Not Calculated'

# Make unknown run types as other
df_facets['FACETS_RUN_TYPE'] = df_facets['FACETS_RUN_TYPE'].fillna('Other')

df_all = df_facets.merge(right=df_cna_portal[['Sample ID', 'CNA Fraction']],
                         how='inner', left_on='SAMPLE_ID', right_on='Sample ID')

df_all = df_all.drop(columns=['Sample ID', 'Unnamed: 0'])
df_all = df_all.rename(columns={'fraction_cna': 'CNA FACETS', 'CNA Fraction': 'CNA Portal'})





# Find s parameter
regex_key = '(s+\d{2,3})'
regex_func = re.compile(regex_key)
list1 = df_all['FOLDER_NAME_FACETS'].apply(lambda str: len(regex_func.findall(str)))
df_all = df_all.assign(HAS_S_PARAMETER=list1)

# Find "norm"
regex_key = '(norm)'
regex_func2 = re.compile(regex_key)
list2 = df_all['FOLDER_NAME_FACETS'].apply(lambda str: len(regex_func2.findall(str)))
df_all = df_all.assign(HAS_NORM=list2)
df_all.head()

# file counts
facets_counts = df_all.groupby(['SAMPLE_ID'])['FOLDER_NAME_FACETS'].count()

# Select cases
# If no S parameter,

# Save table
fname_out = 'facets_fga_all_20190311_with_fga.csv'
df_all.to_csv(fname_out, index=False)

# df_all[df_all.SAMPLE_ID == 'P-0016149-T01-IM6']
tmp = 0
