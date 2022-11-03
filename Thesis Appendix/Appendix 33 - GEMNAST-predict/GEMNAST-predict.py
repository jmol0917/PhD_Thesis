"""
Created by Juan M.
07/2022
"""

import pandas as pd

# read patient data
path = 'C:/'
patient_data = pd.read_csv(path + 'mapped_abundance_general.csv', index_col=0)

# file with fibre degradation data on every AGORA strain and intermediate metabolite synthesis from specific fibre 
# degradation
fibre_degradation = pd.read_csv(path + 'polysaccharide_growth_compilation.csv', index_col=0)
intmet_export = pd.read_csv(path + 'intmet_export_compilation.csv', index_col=0)

# file with information about which mono/oligosaccharides are released from degradation of fibre
ch_metabolic_pathways = pd.read_csv(path + 'ch_metabolic_pathways.csv', index_col=0)

# files with amino acid and vitamin data
amino_acid_compilation = pd.read_csv(path + 'amino_acid_condensed.csv', index_col=0)
vitamin_compilation = pd.read_csv(path + 'cofactor_condensed.csv', index_col=0)

# Based on the fibre and intermediate metabolite of interest, this script compiles strains in the patient of interest that
# are capable of 1) degrading such fibre and 2) synthesising such metabolite from products released from the fibre of 
# interest and reveals amino acids and vitamins required to support the resulting 'network'.

patient_code = ''
fibre = '' # Select fibre from 'polysaccharide_growth_compilation.csv' file. Has to be capitalised.
int_met = ''  # Select one of the following: Acetic acid, Acetaldehyde, Butyrate, Ethanol, Formate, L-Lactic acid, Malic acid, Succinate Propionate, Pyruvate, Succinate
abundance_threshold = 0.0 # Strains with abundances lowe than the specified are not included in final table (main output)

final_dict = {'Model': [], 'Function': [], 'ASVs_represented': [], 'Total_abundance': [], 
              'B1': [], 'B2': [], 'B3': [], 'B5': [], 'B6': [], 'B9': [], 'Kx': [],
              'Glu': [], 'Asp': [], 'Bra': [], 'Ser': [], 'Aro': [], 'His': []
             }

# Here we build a dictionary with tags and models
tag_model_data = {}

for tag, data in patient_data.iterrows():
    if tag not in tag_model_data:
        tag_model_data.update({tag: data['Model']})
        
# Here we are going to identify which mono/oligosaccharides are released from the degradation of the fibre of interest
for substrate, components in ch_metabolic_pathways.iterrows():
    if substrate == fibre:
        simpler_carbohydrates = components[0].split(', ')

# First we identify ASVs mapped to AGORA strains 
# We are going to store the identified AGORA strains in this list
ASVs_in_patient = {}

# We check patient_data file
for patient, reads in patient_data.iteritems():
    # find the patient we are looking for
    if patient == patient_code:
        # find and store AGORA strains within that patient
        for i, read in enumerate(reads):
            if read != 0:
                tag = reads.index[i]
                model = ''
                for entry in tag_model_data:
                    if tag == entry:
                        model = tag_model_data[tag]

                if tag not in ASVs_in_patient and 'Unidentified' not in tag:
                    ASVs_in_patient[tag] = [model, read]


# now that we have a list of the AGORA strains within the patient of interest we can assess their metabolic attributes
# in terms of fibre of interest degradative capabilities and intermediate metabolite degradation
role = {'fibre degrader':[], 'intmet producer': []}

for label, data in fibre_degradation.iteritems():
    if label[2:-3] == fibre:
        for tag in ASVs_in_patient:
            model = ASVs_in_patient[tag][0]
            for strain in data.index:
                if strain == model:
                    if data[strain] == 1:
                        role['fibre degrader'].append([tag, ASVs_in_patient[tag][1]])

for simple_carb in simpler_carbohydrates:            
    for pair, data2 in intmet_synthesis.iteritems():
        if pair == str(simple_carb + ', ' + int_met):
            for tag in ASVs_in_patient:
                model = ASVs_in_patient[tag][0]
                for strain in data2.index:
                    if strain == model:
                        if data2[strain] == 1:
                            if [tag, ASVs_in_patient[tag][1]] not in role['intmet producer']:
                                role['intmet producer'].append([tag, ASVs_in_patient[tag][1]])

# after identifying which ASVs are relevant for the fibre and int met of interest we make a list with them
# regardless of the role attached to them before
models_by_role = {'fibre degrader':[], 'intmet producer': []}
final_strains_of_interest = {}

for asv_role in role:
    role_abundance = 0
    role_strains = []
    for pair in role[asv_role]:
        role_abundance += pair[1]
        role_strains.append(pair)

        if tag_model_data[pair[0]] not in models_by_role[asv_role]:
            models_by_role[asv_role].append(tag_model_data[pair[0]])
        if tag_model_data[pair[0]] not in final_strains_of_interest:
            final_strains_of_interest.update({tag_model_data[pair[0]]: [1 , pair[1]]})
        else:
            final_strains_of_interest[tag_model_data[pair[0]]][0] += 1
            final_strains_of_interest[tag_model_data[pair[0]]][1] += pair[1]
    
    # Total abundance of primary and secondary degraders is printed
    if asv_role == 'fibre degrader':
        print("Total degraders' (primary degraders) abundance:", "{:.4f}".format(role_abundance))
        print("High abundance degraders:")
        for s in role_strains:
            if s[1] > 0.02:
                print(s[0] + ':', s[1])
    if asv_role == 'intmet producer':
        print("Total producers' (secondary degraders) abundance:", "{:.4f}".format(role_abundance))
        print("High abundance producers:")
        for s in role_strains:
            if s[1] > 0.02:
                print(s[0] + ':', s[1])

# Final table assembly strain by strain
# Strains with an abundance lower than the specified in abundance_threshold are not included                 
for strain in final_strains_of_interest:
    if final_strains_of_interest[strain][1] > abundance_threshold:
        final_dict['Model'].append(strain)
        final_dict['ASVs_represented'].append(final_strains_of_interest[strain][0]) 
        final_dict['Total_abundance'].append(final_strains_of_interest[strain][1])

        if strain in models_by_role['intmet producer'] and strain in models_by_role['fibre degrader']:
            final_dict['Function'].append('Both')
        else: 
            if strain in models_by_role['intmet producer']:
                final_dict['Function'].append('Intmet producer')
            if strain in models_by_role['fibre degrader']:
                final_dict['Function'].append('Fibre degrader')
        # here we assign amino acid and vitamin attributes from strain in question
        for index, row in amino_acid_compilation.iterrows():
            if index == strain:
                amino_acids_temp = {'Glu':'', 'Asp':'', 'Bra':'', 'Ser':'', 'Aro':'', 'His':'',} 
                for l, value in enumerate(row):
                    if value == 1:
                        amino_acids_temp[row.index[l][:-4]] = row.index[l][4:]

                for group in amino_acids_temp:
                    final_dict[group].append(amino_acids_temp[group])
        for index2, row2 in vitamin_compilation.iterrows():
            if index2 == strain:
                vitamin_temp = {'B1': '', 'B2': '', 'B3': '', 'B5': '', 'B6': '', 'B9': '', 'Kx': ''}
                for m, value2 in enumerate(row2):
                    if value2 == 1 and 'B12' not in row2.index[m]:
                        vitamin_temp[row2.index[m][:-4]] = row2.index[m][3:]

                for family in vitamin_temp:
                    final_dict[family].append(vitamin_temp[family])
    
fibre_intmet_df = pd.DataFrame.from_dict(final_dict)
fibre_intmet_df.to_csv(path + patient_code + '_' + fibre + '.csv')   
                       
                       
