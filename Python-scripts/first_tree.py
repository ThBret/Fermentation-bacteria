import pandas as pd
import shutil
import os

def get_most_complete_genomes(Acetic_or_Lactic = 'Acetic', data = 'Bacteria log v2.xlsx', path = '/Users/thibaultbret/', Strict = True, Complete = True):
    global most_complete_genomes
    global rows_to_drop
    '''
    Get the most complete genome for each species selected to construct the tree.
    '''
    #make sure the path is the correct one
    if 'Python' in path: path = path.replace('Python','')

    #retrieve data to find the completeness of every genome
    bacteria_log = pd.read_excel(path + data)
    Acetic = bacteria_log[bacteria_log['Lactic/acetic'] == Acetic_or_Lactic]
    most_complete_genomes = Acetic.sort_values('Completeness', ascending = False).drop_duplicates('Unnamed: 0', keep = 'first')[['Unnamed: 0','Assembly Accession']].set_index('Unnamed: 0')
    #if Strict is set to True, remove certain species for the final tree
    if Strict == True:
        to_drop = ['Gluconobacter albidus','Gluconobacter cerevisiae','Gluconobacter japonicus','Gluconobacter oxydans','Neokomagataea thailandica','Neokomagataea tanensis','Neokomagataea anthophila']
        rows_to_drop = most_complete_genomes.loc[to_drop].index
        most_complete_genomes.drop(rows_to_drop, inplace = True)

    #return the assembly accession of the most complete genomes without moving the files
    if Complete != True: return most_complete_genomes

    #get most complete genome files from the Acetic/Lactic directory
    bacteria = [file for file in sorted(os.listdir(path + Acetic_or_Lactic)) if not file.startswith('.')]

    #create new directory containing only the most complete genome per species
    if not os.path.exists(path + Acetic_or_Lactic + '_unique'): os.makedirs(path + Acetic_or_Lactic + '_unique')

    #move the most complete genome files to the new directory
    for bact in bacteria:
        for assemb in most_complete_genomes['Assembly Accession']:
            if bact.startswith(assemb):
                shutil.copyfile(path + Acetic_or_Lactic + '/' + bact, path + '/' + Acetic_or_Lactic + '_unique/' + bact)


if __name__ == "__main__":
    get_most_complete_genomes(data = 'bacteria_log2.xlsx', Complete = False, Strict = False)

