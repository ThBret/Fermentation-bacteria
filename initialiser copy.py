import pandas as pd
from tqdm import tqdm
import os

#prepare the working directory
bacteria_log = pd.DataFrame()
path = '/Users/thibaultbret/Genomes/'
os.chdir(path)
bacteria = [file for file in sorted(os.listdir(path)) if not file.startswith('.')]

for bacterium in bacteria:
    #select the working directory
    path = '/Users/thibaultbret/Genomes/' + bacterium + '/'
    os.chdir(path)
    #fill the bacteria log dataframe using the "data_summary.tsv" file present in the folder downloaded from NCBI
    info_df = pd.read_table('data_summary.tsv', usecols = [0,2,3,4,5,7,8,9,10,11,12,13,14]).rename(columns = {'Organism Scientific Name':'Organism'}).set_index('Organism')
    bacteria_log = bacteria_log.append(info_df)
    #remove unnecessary files
    os.remove('dataset_catalog.json')
    os.remove('assembly_data_report.jsonl')
    #move all the genome sequences to the bacterium-specific folder
    for elem in tqdm(os.listdir(path)):
        if os.path.isdir(elem):
            os.system('mv ' + elem + '/*.fna ' + path)
            if os.path.exists(elem + '/sequence_report.jsonl'): os.remove(elem + '/sequence_report.jsonl')
            os.rmdir(elem)

#Save the full bacteria log data as an Excel sheet
bacteria_log.index = [idx[0] + ' ' + idx[1] for idx in bacteria_log.index.str.split(' ')]
bacteria_log.to_excel('/Users/thibaultbret/bacteria_log00.xlsx')