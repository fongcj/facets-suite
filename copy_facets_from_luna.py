"""
extract_facets_from_luna.py

By Chris Fong - MSKCC 2019


"""
import os
import pandas as pd
from shutil import copyfile


# Data type selection
dtype_suffix = '.seg'

# Read sample list file
path = '/Users/fongc2/Documents/github/MSK/facets-suite/data'
fname = 'sample_ids_aacr_2019.csv'
pathfilename = os.path.join(path, fname)
sample_ids = pd.read_csv(pathfilename, header=0)

# Read directory names, or determine directory list using os.listdir
# Mount example:
path_root = '/Users/fongc2/Desktop/luna_transfer/all/'         # Path for mounted Luna drive
fname_folder_names = 'facet_folder_names_chai.csv'
path_folder_names = ''
pathfilename_folders = os.path.join(path_folder_names, fname_folder_names)
df_folder_names = pd.read_csv(pathfilename_folders, header=0, low_memory=False)

# Define destination folder path
folder_dest = '/Users/fongc2/Desktop/Rdata_subhi_2018_02_27'
if os.path.isdir(folder_dest):
    os.rmdir(folder_dest)
os.makedirs(folder_dest)

# Convert foldername df to list
folder_names = list(df_folder_names['subfolder'])

for i, current_sample_id in enumerate(list(list(sample_ids.iloc[:,0]))):
    matching = [s for s in folder_names if current_sample_id in s]
    if len(matching) == 1:
        # Get list of dir from this folder
        path1 = os.path.join(path_root, matching[0])
        folder_name1 = os.listdir(path1)
        matching1 = [s for s in folder_name1 if 'facets' in s]
        for j, current_facet_dir in enumerate(matching1):
            path2 = os.path.join(path1, current_facet_dir)
            folder_name2 = os.listdir(path2)
            matching2 = [s for s in folder_name2 if dtype_suffix in s]
            # Create folder for each sample
            folder_sample_dest = os.path.join(folder_dest, current_sample_id, current_facet_dir)
            if not os.path.isdir(folder_sample_dest):
                os.makedirs(folder_sample_dest)

            for k, current_rdata in enumerate(matching2):
                #src
                src = os.path.join(path2, current_rdata)
                #dst
                dst = os.path.join(folder_sample_dest, current_rdata)
                # Copy Rdata files over
                print('Copying %s' % current_rdata)
                copyfile(src, dst)
    else:
        # For debugging - will reach this if number of sample folders is not 1
        print('------------------------------------------------------------------------------------------------')
        fdsafas= 0
