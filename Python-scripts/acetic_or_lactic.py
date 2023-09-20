import pandas as pd
import shutil
import os

def write_sort_file(path = os.path.expanduser('~')):
    '''Writes a file ("sort_file.txt") that lists all the genomes and indicates whether each genome is part of the Lactic
    or the Acetic bacteria family.
    '''
    #gather the data 
    bacteria_log = pd.read_excel(os.path.join(path, 'Bacteria log v2.xlsx'))
    sort_dict = bacteria_log.set_index('Assembly Accession').to_dict()['Lactic/acetic']
    sort_dict = {key.split('.')[0]: value for key, value in sort_dict.items()}
    #write the file
    with open(os.path.join(path, 'sort_file.txt'), 'w') as f:
        for k, v in sort_dict.items(): 
            f.write('%s:%s\n' % (k, v))


def reorganiser(path = os.path.expanduser('~'), directory = 'Genomes'):
    '''
    Moves all the ".fna" sequence files from nested directories to the "Genomes" directory.
    '''
    species = [f for f in sorted(os.listdir(os.path.join(path, directory))) if not f.startswith('.')]
    for sp in species:
        genomes = [gen for gen in sorted(os.listdir(os.path.join(path, directory, sp))) if gen.endswith('.fna')]
        for gen in genomes:
            shutil.move(os.path.join(path, directory, sp, gen), os.path.join(path, directory))


def assigner(path = os.path.expanduser('~'), directory = 'Genomes'):
    '''
    Reads the "sort_file.txt" file that assigns every genome to either the Lactic or the Acetic bacteria family and correspondingly
    moves the genome sequence files in the "Acetic" folder or the "Lactic" folder. That way, lactic and acetic bacteria can be treated separately
    for later steps of the analysis.
    '''
    genomes = [file for file in sorted(os.listdir(os.path.join(path, directory))) if not file.startswith('.')]
    if not os.path.exists(os.path.join(path, 'Acetic')): os.makedirs(os.path.join(path, 'Acetic'))
    if not os.path.exists(os.path.join(path, 'Lactic')): os.makedirs(os.path.join(path, 'Lactic'))
    with open(os.path.join(path, 'sort_file.txt')) as f:
        lines = f.readlines()
        for line in lines:
            id, family = line.strip().split(':')
            for gen in genomes:
                if id == gen.split('.')[0]:
                    shutil.copy(os.path.join(path, directory, gen), os.path.join(path, family))

if __name__ == "__main__":
    write_sort_file()
    reorganiser()
    assigner()


