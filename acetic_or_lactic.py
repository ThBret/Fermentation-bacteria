import pandas as pd
import os

def write_sort_file(path = '/Users/thibaultbret/'):
    '''Writes a file ("sort_file.txt") that lists all the genomes and indicates whether each genome is part of the Lactic
    or the Acetic bacteria family.
    '''
    #gather the data 
    bacteria_log = pd.read_excel(path + 'Bacteria log v2.xlsx')
    sort_dict = bacteria_log.set_index('Assembly Accession').to_dict()['Lactic/acetic']
    sort_dict = {key.split('.')[0]: value for key, value in sort_dict.items()}
    #write the file
    with open(path + 'sort_file.txt', 'w') as f:
        for k, v in sort_dict.items(): 
            f.write('%s:%s\n' % (k, v))

def assigner(path = os.getcwd()):
    '''
    Reads the "sort_file.txt" file that assigns every genome to either the Lactic or the Acetic bacteria family and correspondingly
    moves the genome sequence files in the "Acetic" folder or the "Lactic" folder. That way, lactic and acetic bacteria can be treated separately.
    '''
    genomes = [file for file in sorted(os.listdir(path + '/Genomes_filtered')) if not file.startswith('.')]
    if not os.path.exists(path + '/Acetic'): os.makedirs(path + '/Acetic')
    if not os.path.exists(path + '/Lactic'): os.makedirs(path + '/Lactic')
    with open(path + '/sort_file.txt') as f:
        lines = f.readlines()
        for line in lines:
            id, family = line.split(':')
            for genome in genomes:
                if id == genome.split('.')[0]:
                    os.system(f'mv {path}/Genomes_filtered/{genome} {path}/{family}')
