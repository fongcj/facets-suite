"""
facets_folder_struct.py

By Chris Fong - MSKCC 2019


"""
import os
import pandas as pd
import numpy as np
from shutil import copyfile


desired_width = 320
pd.set_option('display.width', desired_width)
np.set_printoptions(linewidth=desired_width)
pd.set_option('display.max_rows', 100)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.expand_frame_repr', False)

# Data type selection
dtype_suffix = '.Rdata'
subfolder_id1 = 'facets'

# Read directory names, or determine directory list using os.listdir
# Mount example:
path_root = '/Users/fongc2/Desktop/luna_transfer/all/'         # Path for mounted Luna drive
# path_root = '/ifs/res/taylorlab/impact_facets/'
fname_folder_names = 'facet_folder_names_chai.csv'
path_folder_names = 'data'
pathfilename_folders = os.path.join(path_folder_names, fname_folder_names)
df_folder_names = pd.read_csv(pathfilename_folders, header=0, low_memory=False)
# folder_names = os.listdir(path_root)
# df_folder_names = pd.DataFrame(folder_names, columns=['subfolder'])
# df_folder_names.to_csv(fname_folder_names, index=False)

# Create data frame containing path and file info, parameters, etc
cols_df = ['FOLDER_NAME_TOP', 'FOLDER_NAME_FACETS', 'FILENAME']
df = pd.DataFrame(columns=cols_df)
df.head()

# Convert foldername df to list
folder_names = list(df_folder_names['subfolder'])
ctr = 0

for i, current_sample_dir in enumerate(df_folder_names.iloc[:, 0]):
    # Get list of dir in sample folder
    path1 = os.path.join(path_root, current_sample_dir)
    try:
        folder_name1 = os.listdir(path1)

        # Find facets folder
        matching1 = [s for s in folder_name1 if subfolder_id1 in s]

        for j, current_facet_dir in enumerate(matching1):
            # Get list of dir in facets folder
            path2 = os.path.join(path1, current_facet_dir)
            try:
                folder_name2 = os.listdir(path2)
                matching2 = [s for s in folder_name2 if dtype_suffix in s]

                for k, current_rdata in enumerate(matching2):
                    # Build rows in dataframe here
                    df.at[ctr] = [current_sample_dir, current_facet_dir, current_rdata]
                    ctr += 1
            except:
                print('No such subdirectory: %s' % path2)
    except:
        print('No such file or directory: %s' % path1)


# Create column for labeling of purity or hisens
is_hisens = df.FILENAME.str.contains('_hisens')
is_purity = df.FILENAME.str.contains('_purity')

df.loc[is_hisens, 'FACETS_RUN_TYPE'] = 'hisens'
df.loc[is_purity, 'FACETS_RUN_TYPE'] = 'purity'


df.to_csv('facets_file_locations_AACR_2019.csv', index=False)

tmp = 0



