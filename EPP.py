import importlib

# Install required packages
packages = ['pandas', 'numpy', 'scipy', 'scikit-posthocs', 'peptides']
for p in packages:
    try:
        importlib.import_module(p)
        print(f"{p} is already installed.")
    except ImportError:
        print(f"{p} is not installed - now installing...")
        try:
            import subprocess
            subprocess.check_call(['pip', 'install', p])
            print(f"{p} has been successfully installed.")
        except Exception as e:
            print(f"Failed to install {p}: {e}")

# Import and abbreviate the packages
import pandas as pd
import numpy as np
import re
from scipy.stats import shapiro, f_oneway, kruskal, levene, ttest_ind, mannwhitneyu
from scikit_posthocs import posthoc_dunn, posthoc_tukey_hsd
from peptides import Peptide


# Writes the core allele dictionary of control and risk HLA alleles
def allele_types(control_alleles, risk_alleles):
    allele_dictionary = {}
    for c in control_alleles:
        allele_dictionary[c] = 'Control'
    for r in risk_alleles:
        allele_dictionary[r] = 'Risk'
    return allele_dictionary


# Formats the IEDB output and creates MixMHC2pred input
def iedb_format(iedb_csv, plist, allele_dictionary):
    data = pd.read_csv(iedb_csv)
    data.rename(columns={'rank': 'adjusted_rank'}, inplace=True)
    filter = (data['adjusted_rank'] <= 20)
    fd = data[filter]
    fd = fd.rename(columns={'seq_num': 'source_protein'})
    values = fd['source_protein'].unique()
    sorted_values = sorted(values)
    fd['source_protein'].replace(sorted_values, plist, inplace=True)
    fd['allele_type'] = fd['allele'].map(allele_dictionary)
    variable_column = fd.pop('allele_type')
    fd.insert(1, 'allele_type', variable_column)
    list = fd.peptide.unique()
    np.asarray(list)
    np.savetxt('HLAII_peptide_output.txt', list, fmt='%s')
    return fd


# Formats the MixMHC2pred output and merges the output
def merge_data(iedb_csv, mix_csv, allele_dictionary):
    with open(mix_csv, "r") as mix_input:
        lines = mix_input.readlines()[19:]
    with open('mix_csv_formatted', "w") as output_file:
        output_file.writelines(lines)
    data2 = pd.read_csv('mix_csv_formatted', delimiter="\t")
    data2['BestAllele_type'] = data2['BestAllele'].map(allele_dictionary)
    best_type_column = data2.pop('BestAllele_type')
    data2.insert(3, 'BestAllele_type', best_type_column)
    data2 = data2.rename(columns={'Peptide': 'peptide'})

    merged_data = pd.merge(iedb_csv, data2, on='peptide')
    merged_data.sort_values('%Rank_best', ascending=True, inplace=True)
    protein_list = merged_data.source_protein.unique()
    DfD = {elem: pd.DataFrame() for elem in protein_list}
    for key in DfD.keys():
        DfD[key] = merged_data[:][merged_data.source_protein == key]
        DfD[key] = DfD[key].sort_values(['start', 'end', 'adjusted_rank'])
    return merged_data, DfD


def pdif(DfD, allele_dictionary):
    def hla_recognise(ad):
        # Create dictionaries for storing the different HLA-II alleles
        dp_strings = {}
        dm_strings = {}
        do_strings = {}
        dq_strings = {}
        dr_strings = {}
        # Create regex expressions for the subtypes
        dp_pattern = re.compile(r'DP')
        dm_pattern = re.compile(r'DM')
        do_pattern = re.compile(r'DO')
        dq_pattern = re.compile(r'DQ')
        dr_pattern = re.compile(r'DR')
        # Execute pattern recognition and add to dictionaries
        for key, value in ad.items():
            if re.search(dp_pattern, key):
                dp_strings[key] = value
            elif re.search(dm_pattern, key):
                dm_strings[key] = value
            elif re.search(do_pattern, key):
                do_strings[key] = value
            elif re.search(dq_pattern, key):
                dq_strings[key] = value
            elif re.search(dr_pattern, key):
                dr_strings[key] = value
        return dp_strings, dm_strings, do_strings, dq_strings, dr_strings
    dp, dm, do, dq, dr = hla_recognise(allele_dictionary)

    subtype_dict = {'dp': dp, 'dm': dm, 'do': do, 'dq': dq, 'dr': dr}
    test_alleles = {subtype: alleles for subtype, alleles in subtype_dict.items() if len(alleles) >= 2}
    pdif_output = pd.DataFrame(columns=['protein'])
    pdif_phs = {}  # Dictionary for post-hocs

    for subtype, alleles in test_alleles.items():
        allele_list = list(alleles.keys())
        if f'{subtype}_test' not in pdif_output.columns:
            pdif_output[f'{subtype}_test'] = None
            pdif_output[f'{subtype}_p'] = None
            pdif_output[f'{subtype}_sig'] = None

        for key in DfD:
            data = DfD[key].loc[DfD[key]['allele'].isin(allele_list), ['adjusted_rank', 'allele']]
            if pdif_output[pdif_output['protein'] == key].empty:
                pdif_output = pdif_output.append({'protein': key}, ignore_index=True)
            temp = []
            violated = False
            for allele in allele_list:
                column_data = data[data['allele'] == allele]['adjusted_rank']
                if len(column_data) >= 3:  # Test how many entries there are in the dataset
                    temp.append(column_data)
                    stat, p = shapiro(column_data)
                    if p >= 0.05:
                        violated = True
                else:
                    temp.append(column_data)
                    violated = True
            if len(allele_list) > 2:
                if not violated:
                    stat, levene_p = levene(*temp)  # Conduct levenes test if data is normal
                    if levene_p < 0.05:
                        violated = True
                if not violated:  # Conduct ANOVA if normal & homogenous variances
                    stat, p = f_oneway(*temp)
                    pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_test'] = 'ANOVA'
                    pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_p'] = p
                    if p <= 0.05:
                        pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_sig'] = 'significant'
                        posthoc = posthoc_tukey_hsd(temp)
                        posthoc = posthoc.rename(columns=dict(zip(posthoc.columns, allele_list)),
                                                 index=dict(zip(posthoc.index, allele_list)))
                        posthoc['total'] = posthoc.apply(lambda row: sum(row <= 0.05), axis=1)
                        pdif_phs[key] = posthoc
                    else:
                        pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_sig'] = 'ns'
                else:  # Kruskal-Wallis if criteria violated
                    stat, p = kruskal(*temp)
                    pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_test'] = 'Kruskal'
                    pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_p'] = p
                    if p <= 0.05:
                        pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_sig'] = 'significant'
                        posthoc = posthoc_dunn(temp, p_adjust='bonferroni')
                        posthoc = posthoc.rename(columns=dict(zip(posthoc.columns, allele_list)),
                                                 index=dict(zip(posthoc.index, allele_list)))
                        pdif_phs[key] = posthoc
                    else:
                        pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_sig'] = 'ns'
            else:
                if not violated:
                    stat, levene_p = levene(*temp)  # Conduct levenes test if data is normal
                    if levene_p < 0.05:
                        violated = True
                if not violated:  # Conduct T-test if normality satisfied
                    stat, p = ttest_ind(*temp)  # T-test if criteria satisfied
                    pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_test'] = 't_test'
                    pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_p'] = p
                    if p <= 0.05:
                        pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_sig'] = 'significant'
                    else:
                        pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_sig'] = 'ns'
                else:
                    stat, p = mannwhitneyu(*temp)  # MannwhitneyU if criteria violated
                    pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_test'] = 'mannwhitneyu'
                    pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_p'] = p
                    if p <= 0.05:
                        pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_sig'] = 'significant'
                    else:
                        pdif_output.loc[pdif_output['protein'] == key, f'{subtype}_sig'] = 'ns'

    return pdif_output, pdif_phs


def adif(merged_data):
    allele_list = merged_data.allele.unique()
    AD = {elem: pd.DataFrame() for elem in allele_list}
    for key in AD.keys():
        AD[key] = merged_data[:][merged_data.allele == key]

    temp = []
    adif_output = pd.DataFrame(columns=['Allele', 'test', 'p_value', 'sig'])
    adif_phs = {}

    for key in AD:
        data = AD[key].loc[:, ['adjusted_rank', 'source_protein']]
        adif_output = adif_output.append({'Allele': key}, ignore_index=True)
        proteins = AD[key].source_protein.unique()
        violated = False

        for protein in proteins:
            column_data = data[data['source_protein'] == protein]['adjusted_rank']
            temp.append(column_data)
            if len(column_data) < 3:
                violated = True
            else:
                stat, p = shapiro(column_data)
                if p >= 0.05:
                    violated = True
        if not violated:
            stat, levene_p = levene(*temp)
            if levene_p < 0.05:
                violated = True
        if not violated and len(proteins) > 2:
            stat, p = f_oneway(*temp)  # ANOVA if criteria satisfied
            adif_output.loc[adif_output['Allele'] == key, 'test'] = 'ANOVA'
            if p <= 0.05:
                posthoc = posthoc_tukey_hsd(temp)
                posthoc = posthoc.rename(columns=dict(zip(posthoc.columns, proteins)),
                                         index=dict(zip(posthoc.index, proteins)))
                posthoc['total'] = posthoc.apply(lambda row: sum(row <= 0.05), axis=1)
                adif_phs[key] = posthoc
        if not violated and len(proteins) <= 2:
            stat, p = ttest_ind(*temp)  # T-test if criteria satisfied
            adif_output.loc[adif_output['Allele'] == key, 'test'] = 't_test'
        if violated and len(proteins) > 2:
            stat, p = kruskal(*temp)  # Kruskal-Wallis if criteria violated
            adif_output.loc[adif_output['Allele'] == key, 'test'] = 'Kruskal'
            if p <= 0.05:
                posthoc = posthoc_dunn(temp, p_adjust='fdr_bh')
                posthoc = posthoc.rename(columns=dict(zip(posthoc.columns, proteins)),
                                         index=dict(zip(posthoc.index, proteins)))
                posthoc['total'] = posthoc.apply(lambda row: sum(row <= 0.05), axis=1)
                adif_phs[key] = posthoc
        if violated and len(proteins) <= 2:
            stat, p = mannwhitneyu(*temp)
            adif_output.loc[adif_output['Allele'] == key, 'test'] = 'mannwhitneyu'
        temp = []

        adif_output.loc[adif_output['Allele'] == key, 'p_value'] = p
        if p <= 0.05:
            adif_output.loc[adif_output['Allele'] == key, 'sig'] = 'sig'
        else:
            adif_output.loc[adif_output['Allele'] == key, 'sig'] = 'ns'
        adif_output.sort_values('p_value', ascending=True, inplace=True)
        for protein in proteins:
            temp_mean = data[data['source_protein'] == protein]['adjusted_rank'].mean()
            protein_col = protein + '_mean'
            adif_output.loc[adif_output['Allele'] == key, protein_col] = temp_mean

    return adif_output, adif_phs


def identify_regions(merged_data, target_proteins, allele_dictionary, length=15, abv_proportion=0.7, peptide_pH=7.4):
    targets = [target_proteins]
    length -= 1
    regions = {elem: pd.DataFrame() for elem in targets}
    peptide = {elem: pd.DataFrame() for elem in targets}

    class hla_allele:
        def __init__(self, a, b, c):
            self.subtype = a
            self.name = b
            self.type = c

    def hla_recognise(ad):
        allele_list = []  # Initialize as a list
        # Create regex expressions for the subtypes
        dp_pattern = re.compile(r'DP')
        dm_pattern = re.compile(r'DM')
        do_pattern = re.compile(r'DO')
        dq_pattern = re.compile(r'DQ')
        dr_pattern = re.compile(r'DR')
        # Execute pattern recognition and add to list
        for key, value in ad.items():
            if re.search(dp_pattern, key):
                subtype = 'dp'
            elif re.search(dm_pattern, key):
                subtype = 'dm'
            elif re.search(do_pattern, key):
                subtype = 'do'
            elif re.search(dq_pattern, key):
                subtype = 'dq'
            elif re.search(dr_pattern, key):
                subtype = 'dr'
            temp = hla_allele(subtype, key, value)
            allele_list.append(temp)
        return allele_list

    allele_list = hla_recognise(allele_dictionary)
    for alleles in allele_list:
        for target in targets:
            # 1 Identify protein regions
            region_output = pd.DataFrame(columns=['RegionID', 'RegionStart', 'RegionEnd', 'allele', 'allele_type'])
            if alleles.type == 'Risk':
                indicator = 'r'
            else:
                indicator = 'c'
            data = merged_data[merged_data['allele'] == alleles.name]

            RegionID = 0
            RegionStart = data['start'].min()
            RegionEnd = RegionStart + length
            while True:
                subset = data.loc[(data['start'] >= RegionStart) & (data['start'] <= RegionEnd), 'start']
                if not subset.empty:
                    max_start = max(subset)
                    new_RegionEnd = max_start + length
                    if new_RegionEnd <= RegionEnd:
                        RegionID += 1
                        region_output = region_output.append(
                            {'RegionID': f'{target}_{indicator}{RegionID}', 'RegionStart': RegionStart,
                             'RegionEnd': RegionEnd, 'allele': alleles.name, 'allele_type': alleles.type},
                            ignore_index=True)
                        RegionStart = data.loc[data['start'] > RegionEnd, 'start'].min()
                        RegionEnd = RegionStart + length
                        continue
                    else:
                        RegionEnd = new_RegionEnd
                else:
                    break
            target_dict = regions.setdefault(target, {})
            subtype_list = target_dict.setdefault(alleles.subtype, [])
            subtype_list.append(region_output)
    return regions


"""
 # 2 Identify candidate peptides using protein ABVs
        peptide_output = pd.DataFrame()
        db_peptides = DfDRB[target]
        for index, row in db_peptides.iterrows():
            if row['allele_type'] == 'Risk':
                risk_peptide = row['peptide']
                control_row = db_peptides[
                    (db_peptides['allele_type'] == 'Control') & (db_peptides['peptide'] == risk_peptide)]
                if control_row.empty:
                    peptide_output = peptide_output.append(row)
                else:
                    risk_adjusted_rank = row['adjusted_rank']
                    control_adjusted_rank = control_row['adjusted_rank'].iloc[0]
                    if risk_adjusted_rank <= abv_proportion * control_adjusted_rank:
                        peptide_output = peptide_output.append(row)
        file_name = 'DRB_' + str(target) + '_peptides.txt'
        peptide_output['peptide'].to_csv(f'Peptides/DRB/{file_name}', index=False, header=True)
        # 3 Predict the peptide properties using peptides.py
        pep_df = pd.DataFrame(columns=['peptide', 'charge', 'hydrophobicity', 'instability',
                                       'isoelectric_point', 'molecular_weight', 'mz'])
        for index, row in peptide_output.iterrows():
            peptide = Peptide(row['peptide'])
            pep_df.loc[index] = [row['peptide'], peptide.charge(pH=peptide_pH, pKscale="EMBOSS"),
                                 peptide.hydrophobicity(scale="Aboderin"), round(peptide.instability_index(), 2),
                                 peptide.isoelectric_point(pKscale="EMBOSS"), peptide.molecular_weight(), peptide.mz()]
        peptide[target] = pep_df
        regions[target] = peptide_output
        regions[target].sort_values('adjusted_rank', ascending=True, inplace=True)

"""