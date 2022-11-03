import pandas as pd
import os
import json
import itertools
from os import listdir
from os.path import isfile, join

# file with fibre degradation data on every AGORA strain and intermediate metabolite synthesis from specific fibre 
# degradation
path = 'C:/'

# carbohydrate degradation and intermediate generation files
poly_degradation = pd.read_csv(path + 'polysaccharide_growth_compilation.csv', index_col=0)
proteoglycan_deg = pd.read_csv(path + 'proteoglycan_growth_compilation.csv', index_col=0)
monosaccharide_deg = pd.read_csv(path + 'monosaccharide_growth_compilation.csv', index_col=0)
oligosaccharide_deg = pd.read_csv(path + 'oligosaccharide_growth_compilation.csv', index_col=0)
intmet_export = pd.read_csv(path + 'intmet_export_internisable.csv', index_col=0)
intmet_deg = pd.read_csv(path + 'intmet_growth_compilation.csv', index_col=0)

# file with information about which mono/oligosaccharides are released from degradation of fibre
ch_metabolic_pathways = pd.read_csv(path + 'ch_metabolic_pathways.csv', index_col=0)

# files with amino acid and vitamin data
amino_acid_compilation = pd.read_csv(path + 'amino_acid_condensed.csv', index_col=0)
vitamin_compilation = pd.read_csv(path + 'cofactor_condensed.csv', index_col=0)

# network file
# change name for R or B networks
network_file = pd.read_csv(path + 'B_network_with_models.csv', index_col=0)

    
strain_categories = {
    'Mucin': ['Core 2', 'Core 3', 'Core 4', 'Core 5', 'Core 6', 'Core 7', 'Core 8', 
             'Sialyl-T Antigen', 'Sialyl-Tn Antigen', 'T-Antigen (Core 1)', 'Tn Antigen',
             'GlcNac-Alpha-1,4-Core 1', 'GlcNac-Alpha-1,4-Core 2'
             ],
    'Heparan sulphate': ['Heparan Sulfate Proteoglycan', 'Product 1 HEPARL1_e', 'Product AS3TASE_HS2', 
                         'Product GAM2STASE_HS3', 'Product GLCNACASE_HS2', 'Product GLCAASE_HS', 'Product GLCNACASE_HS3',
                         'Product IDOURASE_HS3', 'Product 2 HEPARL1_e', 'Product IDOURASE_HS1', 'Product AS6TASE_HS1',
                         'Product AS3TASE_HS1', 'Product GAM2STASE_HS2', 'Product IS2TASE_HS1', 'Product IDOURASE_HS2',
                         'Product AS6TASE_HS2'],
    'Cartilage': ['Chondroitin Sulfate A Proteoglycan', 'Dimer GalNAc4S-GlcA', 'Chondroitin Sulfate B Proteoglycan',
                  'Dimer GalNAc4S-IdoA2S', 'Chondroitin Sulfate C Proteoglycan', 'Dimer GalNAc6S-GlcA',
                  'Disialyl-T Antigen', 'Hyaluronic acid'],
    'Glycogen': ['Glycogen, Structure 2', 'Glycogen, Structure 4', 'Glycogen, Structure 5'],
    'Starch': ['Starch', 'Starch, Structure 1', 'Amylopectin', 'Amylose'],
    'Inulin': ['Inulin', 'Kestopentaose', 'Kestotetraose'],
    'Other Fibre': ['Arabinan', 'Larch arabinogalactan', 'Arabinoxylan', 'Beta-Glucans', 'Cellulose', 'Dextran', 
                    'Carob galactomannan', 'Homogalacturonan', 'Inulin', 'Levan', 'Lichenin', 'Laminarin', 'Alpha-mannan, yeast', 'Pectin',
                    'Pectic galactan (potato)', 'Pullulan', 'Potato rhamnogalacturonan I', 'Wine rhamnogalacturonan II',
                    'Xylan','Xyluglucan'],
    'Oligo': ['Arabinotriose', 'Cellobiose', 'Dextrin', 'N,N-diacetylchitobiose', 'D-Maltose', 
              'Kestopentaose', 'Kestotetraose', 'Maltohexaose', 'Mannotriose (beta-1,4)', 'Melibiose', 
              'Raffinose', 'Stachyose', 'Sucrose', 'Trehalose', 'Alpha-Lactose', 'Maltotriose'],
    'Mono': ['N-Acetylgalactosamine', 'N-Acetyl-D-glucosamine', 'N-Acetylneuraminic acid', 'D-Arabinose', 
             'L-Arabinose', 'Deoxyribose', 'D-Fructose', 'L-Fucose', 'Glucosamine', 'D-Glucose', 'L-lyxose',
             'D-Mannose', 'D-Ribose', 'L-Rhamnose', 'Salicin', 'D-Xylose', 'D-Galactose', 'D-Galacturonate',
             'D-Glucuronic acid'],
    'Int_met': ['Acetic acid', 'Acetaldehyde', 'Ethanol', 'Formic acid', 'L-Lactic acid', 'Malic Acid' 
                        'Propionate', 'Pyruvic acid', 'Butyrate', 'Succinate', 'Fumarate']
}

extracellular_oligos = ['Kestopentaose', 'Kestotetraose', 'Arabinotriose']


final_dict = {'Tag':[], 'Strain': [], 'Proteoglycans': [], 'Polysaccharides': [], 'CH_generated': [], 
              'Oligosaccharides': [], 'Monosaccharides': [], 'IM_generated': [], 'Int_met': [], 
              'B1': [], 'B2': [], 'B3': [], 'B5': [], 'B6': [], 'B9': [], 'Kx': [],
              'Glu': [], 'Asp': [], 'Bra': [], 'Ser': [], 'Aro': [], 'His': []
             }

network_members = {}

for interaction, details in network_file.iterrows():
    source = details['source_models']
    source_tag = interaction
    target = details['target_models']
    target_tag = details['target']
    
    if source_tag not in network_members and 'Unidentified' not in source:
        network_members.update({source_tag:source})
    if target_tag not in network_members and 'Unidentified' not in target:
        network_members.update({target_tag:target})
# print(network_members)

for member in network_members:
    tag = member
    strain = network_members[member]
    
    # categories to be filled
    proteoglycans = []
    polysaccharides = []
    generated_ch = []
    generated_intmet = []
    oligosaccharides, monosaccharides, intermediate_metabolites = [], [], []
    
    # define which polysaccharides the strain degrades and which oligo/monosaccharides are released
    for strain_poly, fibres in poly_degradation.iterrows():
        if strain == strain_poly:
            for i, fibre in enumerate(fibres):
                if fibre == 1: 
                    if len(fibres.index[i]) > 2:
                        fibre_name = fibres.index[i][2:-3]
                    if fibre_name in strain_categories['Other Fibre']:
                        polysaccharides.append(fibre_name)
                    if fibre_name in strain_categories['Starch'] and 'Starch' not in polysaccharides:
                        polysaccharides.append('Starch')
                    if fibre_name in strain_categories['Glycogen'] and 'Glycogen' not in polysaccharides:
                        polysaccharides.append('Glycogen')
                    for molecule, products in ch_metabolic_pathways.iterrows():
                        if molecule == fibre_name:
                            for product in products[0].split(', '):
                                if product in strain_categories['Mono'] or product in strain_categories['Oligo']:
                                    if product not in generated_ch:
                                        generated_ch.append(product)
    
    # define which proteoglycans the strain degrades and which oligo/monosaccharides are released
    for strain_prot, glycans in proteoglycan_deg.iterrows():
        if strain == strain_prot:
            for i, glycan in enumerate(glycans):
                if glycan == 1:
                    if len(glycans.index[i]) > 2:
                        glycan_name = glycans.index[i][2:-3]
                    if glycan_name in strain_categories['Mucin'] and 'Mucin' not in proteoglycans:
                        proteoglycans.append('Mucin')
                    if glycan_name in strain_categories['Heparan sulphate'] and 'Heparan sulphate' not in proteoglycans:
                        proteoglycans.append('Heparan sulphate')
                    if glycan_name in strain_categories['Cartilage'] and 'Cartilage' not in proteoglycans:
                        proteoglycans.append('Cartilage')
                    if glycan_name in strain_categories['Antigens'] and 'Antigens' not in proteoglycans:
                        proteoglycans.append('Antigens')
                    for molecule, products in ch_metabolic_pathways.iterrows():
                        if molecule == glycan_name:
                            for product in products[0].split(', '):
                                if product in strain_categories['Mono'] or product in strain_categories['Oligo']:
                                    if product not in generated_ch:
                                        generated_ch.append(product)
    
    # define which oligosaccharides the strain degrades
    for strain_oligo, oligos in oligosaccharide_deg.iterrows():
        if strain == strain_oligo:
            for i, oligo in enumerate(oligos):
                if oligo == 1:
                    oligosaccharides.append(oligos.index[i])
                    if oligos.index[i] in extracellular_oligos:
                        for molecule, products in ch_metabolic_pathways.iterrows():
                            if molecule == oligos.index[i]:
                                for product in products[0].split(', '):
                                    if product in strain_categories['Mono'] and product not in generated_ch:
                                        generated_ch.append(product)                 
                    
    # define which monosaccharides the strain degrades
    for strain_mono, monos in monosaccharide_deg.iterrows():
        if strain == strain_mono:
            for i, mono, in enumerate(monos):
                if mono == 1:
                    monosaccharides.append(monos.index[i])
                    
    # define which int_met are synthesised from mono/oligosaccharides
    for strain_intmet_syn, pairs in intmet_export.iterrows():
        if strain == strain_intmet_syn:
            for i, pair in enumerate(pairs):
                if pair == 1:
                    real_pair = pairs.index[i].split(', ')
                    if real_pair[0] in monosaccharides and real_pair[1] in strain_categories['Int_met']:
                        if real_pair[1] not in generated_intmet:
                            generated_intmet.append(real_pair[1])
                    if real_pair[0] in oligosaccharides and real_pair[1] in strain_categories['Int_met']:
                        if real_pair[0] not in extracellular_oligos:
                            if real_pair[1] not in generated_intmet:
                                generated_intmet.append(real_pair[1])
    
    # define which int_met can be utilised for growth
    for strain_intmet_deg, intmets in intmet_deg.iterrows():
        if strain == strain_intmet_deg:
            for i, intmet in enumerate(intmets):
                if intmet == 1 and intmets.index[i] in strain_categories['Int_met']:
                    intermediate_metabolites.append(intmets.index[i])
                    
    # we append all the characteristics from the current strain into the final dictionary's corresponding categories
    final_dict['Tag'].append(tag)
    final_dict['Strain'].append(strain)
    final_dict['Proteoglycans'].append(str(proteoglycans)[1:-1].replace("'", ""))
    final_dict['Polysaccharides'].append(str(polysaccharides)[1:-1].replace("'", ""))
    final_dict['CH_generated'].append(str(generated_ch)[1:-1].replace("'", ""))
    final_dict['IM_generated'].append(str(generated_intmet)[1:-1].replace("'", ""))
    final_dict['Oligosaccharides'].append(str(oligosaccharides)[1:-1].replace("'", ""))
    final_dict['Monosaccharides'].append(str(monosaccharides)[1:-1].replace("'", ""))
    final_dict['Int_met'].append(str(intermediate_metabolites)[1:-1].replace("'", ""))  
    
#     here we assign amino acid and vitamin requirements to the strain in question
    for index, row in amino_acid_compilation.iterrows():
        if index == strain:
            amino_acids_temp = {'Glu':'', 'Asp':'', 'Bra':'', 'Ser':'', 'Aro':'', 'His':'',} 
            for l, value in enumerate(row):
                if value == 1:
                    amino_acids_temp[row.index[l][:-4]] = row.index[l][4:]
#             print(strain, amino_acids_temp)
            for group in amino_acids_temp:
                final_dict[group].append(amino_acids_temp[group])
    for index2, row2 in vitamin_compilation.iterrows():
        if index2 == strain:
            vitamin_temp = {'B1': '', 'B2': '', 'B3': '', 'B5': '', 'B6': '', 'B9': '', 'Kx': ''}
            for m, value2 in enumerate(row2):
                if value2 == 1 and 'B12' not in row2.index[m]:
                    vitamin_temp[row2.index[m][:-4]] = row2.index[m][3:]
#             print(strain, vitamin_temp)
            for family in vitamin_temp:
                final_dict[family].append(vitamin_temp[family])
                        
capabilities_network = pd.DataFrame.from_dict(final_dict)

# Uncomment to generate a file at the individual strain level 
## Change name to B or R based on the network assessed
# capabilities_network.to_csv(path + 'B_network_capabilities.csv') 

# this is the dictionary where we store the data derived from our 'degrader' analysis where we identify which nutrients
# individual strains can provided in the context of the explored network
degrader_perspective = {
    'Tag': [], 'Strain': [], 'Shared_AAs': [], 'Shared_Vits': [], 'Shared_CHs': [], 'Shared_IMs': []
}


# the following dictionaries count and store the number of times a nutrient appears in our network and numbers are reported
# at the end
# some energy sources can appear in both_degraded and shared_products
shared_aminos = [
    'Asp', 'Glu', 'Aro', 'Ser', 'Bra', 'His'
]
shared_vits = [
    'B1', 'B2', 'B3', 'B5', 'B6', 'B9', 'Kx'
]

# here we store strains that have been assessed as degraders so we don't have any duplicates at the end as individual
# strains are likely to appear more than one time (interact with more than one strain) in a network file
assessed_degraders = []

shared_goods = {}

every_pair_back_and_forth = []

for interaction, details in network_file.iterrows():
    every_pair_back_and_forth.append([interaction, details['target'], details['edge_width']])
    every_pair_back_and_forth.append([details['target'], interaction, details['edge_width']])

for pair in every_pair_back_and_forth:
    source_tag = pair[0]
    focus_capabilities = []
    generated_chs = []
    generated_int = []
    essential_vits = []
    vits_not_provided_by_network = []
    vits_for_the_network = []
    optional_vits = []
    essential_aas = []
    aas_not_provided_by_network = []
    aas_for_the_network = []
    optional_aas = []
    if 'Unidentified' not in source_tag and source_tag not in assessed_degraders:
        assessed_degraders.append(source_tag)
        strain = ''
        for member, capabilities in capabilities_network.iterrows():
            if capabilities['Tag'] == source_tag:
                strain = capabilities['Strain']
                focus_capabilities = capabilities
                for i, capability in enumerate(capabilities):
                    if capability == 'req' and capabilities.index[i] in shared_vits:
                        essential_vits.append(capabilities.index[i])
                    if capability == 'syn' and capabilities.index[i] in shared_vits:
                        optional_vits.append(capabilities.index[i])
                    if capability == 'req' and capabilities.index[i] in shared_aminos:
                        essential_aas.append(capabilities.index[i])
                    if capability == 'syn' and capabilities.index[i] in shared_aminos:
                        optional_aas.append(capabilities.index[i])
                    if capability == 'sec' and capabilities.index[i] in shared_aminos:
                        optional_aas.append(capabilities.index[i])
                vits_not_provided_by_network = essential_vits
                aas_not_provided_by_network = essential_aas
#         print(source_tag, optional_vits, essential_vits)
#         print(source_tag, optional_aas, essential_aas)
        
        for pair2 in every_pair_back_and_forth:
            if pair2[0] == source_tag and 'Unidentified' not in pair2[1] :
                partner_tag = pair2[1]
                partner_capabilities = []
                partner_essential_vits = []
                partner_optional_vits = []
                partner_essential_aas = []
                partner_optional_aas = []
                for member2, capabilities2 in capabilities_network.iterrows():
                    if capabilities2['Tag'] == partner_tag:
                        partner_capabilities = capabilities2
                        for i, capability in enumerate(capabilities2):
                                if capability == 'req' and capabilities2.index[i] in shared_vits:
                                    partner_essential_vits.append(capabilities2.index[i])
                                if capability == 'syn' and capabilities2.index[i] in shared_vits:
                                    partner_optional_vits.append(capabilities2.index[i])
                                if capability == 'req' and capabilities2.index[i] in shared_aminos:
                                    partner_essential_aas.append(capabilities2.index[i])
                                if capability == 'syn' and capabilities2.index[i] in shared_aminos:
                                    partner_optional_aas.append(capabilities2.index[i])
                                if capability == 'sec' and capabilities2.index[i] in shared_aminos:
                                    partner_optional_aas.append(capabilities2.index[i])
#                 print(partner_tag, partner_optional_vits, partner_essential_vits)
#                 print(partner_tag, partner_optional_aas, partner_essential_aas)

                if pair2[2] > 0:
                    mono_provided = set(focus_capabilities['CH_generated'].split(', ')).intersection(partner_capabilities['Monosaccharides'].split(', '))
                    for mono in mono_provided:
                        if len(mono) > 0:
                            if mono not in shared_goods:
                                shared_goods.update({mono: 1})
                            else:
                                shared_goods[mono] += 1
                        if mono not in generated_chs:
                            generated_chs.append(mono)
                    
                    oligo_provided = set(focus_capabilities['CH_generated'].split(', ')).intersection(partner_capabilities['Oligosaccharides'].split(', '))
                    for oligo in oligo_provided:
                        if len(oligo) > 0:
                            if oligo not in shared_goods:
                                shared_goods.update({oligo: 1})
                            else:
                                shared_goods[oligo] += 1
                        if oligo not in generated_chs:
                            generated_chs.append(oligo)

                    int_provided = set(focus_capabilities['IM_generated'].split(', ')).intersection(partner_capabilities['Int_met'].split(', '))
                    for im in int_provided:
                        if len(im) > 0:
                            if im not in shared_goods:
                                shared_goods.update({im: 1})
                            else:
                                shared_goods[im] += 1
                        if im not in generated_int:
                            generated_int.append(im)
                    
                    vits_for_partner = set(optional_vits).intersection(partner_essential_vits)
                    for vit in vits_for_partner:
                        if len(vit) > 0:
                            if vit not in shared_goods:
                                shared_goods.update({vit: 1})
                            else:
                                shared_goods[vit] += 1
                        if vit not in vits_for_the_network:
                            vits_for_the_network.append(vit)
                    
                    vits_to_self = set(essential_vits).intersection(partner_optional_vits)
                    for vit in vits_to_self:
                        if vit in vits_not_provided_by_network:
                            vits_not_provided_by_network.remove(vit)

                    aas_for_partner = set(optional_aas).intersection(partner_essential_aas)
                    for aa in aas_for_partner:
                        if len(aa) > 0:
                            if aa not in shared_goods:
                                shared_goods.update({aa: 1})
                            else:
                                shared_goods[aa] += 1
                        if aa not in aas_for_the_network:
                            aas_for_the_network.append(aa)
                            
                    aas_to_self = set(essential_aas).intersection(partner_optional_aas)
                    for aa in aas_to_self:
                        if aa in aas_not_provided_by_network:
                            aas_not_provided_by_network.remove(aa)
                    
        degrader_perspective['Tag'].append(source_tag)
        degrader_perspective['Strain'].append(strain)
        degrader_perspective['Shared_AAs'].append(str(sorted(aas_for_the_network)).replace("'", "")[1:-1])
        degrader_perspective['Shared_Vits'].append(str(sorted(vits_for_the_network)).replace("'", "")[1:-1])
        degrader_perspective['Shared_CHs'].append(str(sorted(generated_chs)).replace("'", "")[1:-1])
        degrader_perspective['Shared_IMs'].append(str(sorted(generated_int)).replace("'", "")[1:-1])

with open(path + "shared_goods.txt", 'w') as f: 
    for key, value in shared_goods.items(): 
        f.write('%s:%s\n' % (key, value))

degrader_perspective_df = pd.DataFrame.from_dict(degrader_perspective)
degrader_perspective_df.to_csv(path + 'B_network_public_goods.csv')
