import json
import sys
import pandas as pd
from argparse import ArgumentParser


def get_arguments():
    parser = ArgumentParser(description='Parse mykrobe predict JSON files')

    # job submission options
    parser.add_argument('--jsons', required=True, nargs='+', help='JSON files output from mykrobe predict')
    parser.add_argument('--prefix', required=True, help='prefix for output files')

    return parser.parse_args()


def set_amr_NA(amr_out_dict, drug_linked_mutations):
    for drug in drug_linked_mutations:
        amr_out_dict[drug] = ['NA']
        for gene in drug_linked_mutations[drug]:
            amr_out_dict[gene] = ['NA']
    amr_out_dict['num QRDR'] = ['NA']
    return amr_out_dict


def extract_amr_info(genome_data, genome_name, spp_call):
    # make a dict of all drugs linked with all of their possible mutations
    drug_linked_mutations = {'azithromycin': ['acrB_R717L', 'acrB_R717Q', 'mphA', 'ermB', 'ereA'],
                             'ampicillin': ['blaTEM-1D', 'blaOXA-7'], 'ceftriaxone': ['blaCTX-M-15', 'AmpC1', 'blaOXA-134', 'blaSHV-12'],
                             'ciprofloxacin': ['parC_S80I', 'parC_S80R', 'parC_E84G', 'parC_E84K', 'gyrA_S83F',
                                               'gyrA_S83Y', 'gyrA_D87G', 'gyrA_D87N', 'gyrA_D87V', 'gyrA_D87Y',
                                               'gyrB_S464F', 'gyrB_S464Y', 'qnrS1', 'qnrB1', 'qnrD1'],
                             'chloramphenicol': ['catA1', 'cmlA1'], 'sulfonamides': ['sul1', 'sul2'],
                             'trimethoprim': ['dfrA1', 'dfrA5', 'dfrA7', 'dfrA14', 'dfrA15', 'dfrA17', 'dfrA18'],
                             'tetracycline': ['tetA', 'tetB', 'tetC', 'tetD'],
                             'IncFIAHI1': ['IncFIAHI1'], 'IncHI1A': ['IncHI1A'], 'IncHI1BR27': ['IncHI1BR27'],
                             'IncY': ['IncY'], 'IncX3':['IncX3'], 'IncHI2A': ['IncHI2A'], 'IncI1':['IncI1'],
                             'IncL_M': ['IncL_M'], 'IncFIB_pHCM2': ['IncFIB_pHCM2'], 'IncFIB_K': ['IncFIB_K'],
                             'IncN': ['IncN'], 'z66': ['z66'], 'pST': ['NA_C19241A']}
    # make a dict to convert the beta lactamase enzyme names from mykrobe into their more correct genetic names
    bla_enzymes = {'TEM1D': 'blaTEM-1D', 'CTXM15': 'blaCTX-M-15', 'OXA7': 'blaOXA-7', 'OXA134': 'blaOXA-134', 'SHV12': 'blaSHV-12'}
    # list of all possible plasmid reps so we can make these 0/1 rather than R/S
    plasmid_reps = ['IncFIAHI1', 'IncHI1A', 'IncHI1BR27', 'IncY', 'z66', 'pST', 'IncX3', 'IncHI2A', 'IncI1', 'IncL_M','IncFIB_pHCM2','IncFIB_K','IncN']

    # set up empty dict for final output
    amr_out_dict = {}
    res_data = None
    # create the pandas dataframe for this genome
    amr_out_dict['genome'] = [genome_name]
    # initalise the number of QRDR gyr/par mutations
    num_qrdr_calls = 0
    # intalise sulfamethoxazole variables
    sulfamethoxazole = 0
    sulfamethoxazole_determinants = []
    # can only extract qrdr info if sample is sonnei, otherwise set 'unknown' for all calls
    if spp_call == "Salmonella_Typhi":
        try:
            res_data = genome_data["susceptibility"]
        # if nothing has been called, just set everything to NA
        except KeyError:
            amr_out_dict = set_amr_NA(amr_out_dict, drug_linked_mutations)
    # if the species hasn't been called, then set everything to NA
    else:
        amr_out_dict = set_amr_NA(amr_out_dict, drug_linked_mutations)
    # if no resistance class have been made, then also set everything to NA
    if not res_data:
        amr_out_dict = set_amr_NA(amr_out_dict, drug_linked_mutations)

    else:
        # loop through each drug
        for drug in list(res_data.keys()):
            possible_genes = drug_linked_mutations[drug]
            if res_data[drug]['predict'] == 'R' or res_data[drug]['predict'] == 'r':
                orig_res_calls = res_data[drug]['called_by']
                new_res_calls = []
                # fix any of the gyrA/parC/acrB calls to remove the sequence
                for call in orig_res_calls:
                    call = call.split('-')[0]
                    # if it's one of the bla enzymes, update the name to the genetically correct versio
                    if call in bla_enzymes:
                        new_res_calls.append(bla_enzymes[call])
                    else:
                        new_res_calls.append(call)
                # loop through each call, and record 0/1 depending on whether its present
                for gene in possible_genes:
                    if gene in new_res_calls:
                        amr_out_dict[gene] = [1]
                    else:
                        amr_out_dict[gene] = [0]
                # now added R/S for the drug class column, unless it is cipro in which case we need to do other checks
                if drug == 'ciprofloxacin':
                    # rules for ciprofloxacin are:
                    #1 or 2 QRDR, no qnr -> I
                    #qnr, no QRDR -> I
                    #3 QRDR -> R
                    #>1 QRDR + qnr -> R
                    # count the number of QRDR mutations
                    qnr_present = False
                    for call in new_res_calls:
                        if call.startswith('gyr') or call.startswith('parC'):
                            num_qrdr_calls += 1
                        elif call.startswith('qnr'):
                            qnr_present = True
                    if (num_qrdr_calls == 3) or (num_qrdr_calls >= 1 and qnr_present):
                        res_string = 'R: '
                    elif (not qnr_present and num_qrdr_calls >= 1) or (qnr_present and num_qrdr_calls == 0) or (num_qrdr_calls >=1 and num_qrdr_calls < 3):
                        res_string = 'I: '
                    else:
                        print('THIS SHOULDNT HAPPEN')
                        res_string = 'S: '
                    amr_out_dict['ciprofloxacin'] = [res_string + ';'.join(new_res_calls)]
                # don't do this for the plasmid reps, only for the actual drugs
                elif drug not in plasmid_reps:
                    amr_out_dict[drug] = ['R: ' + ';'.join(new_res_calls)]
                # if we have both a dfr gene and a sul gene, we want to make a sulfamethoxazole column
                # so sulfamethoxazole has to equal 2 for this to become R
                if drug == "trimethoprim":
                    sulfamethoxazole += 1
                    sulfamethoxazole_determinants += new_res_calls
                if drug == "sulfonamides":
                    sulfamethoxazole += 1
                    sulfamethoxazole_determinants += new_res_calls
            # otherwise this genome is susceptible so record it as such
            else:
                amr_out_dict[drug] = ['S']
                for gene in possible_genes:
                    amr_out_dict[gene] = [0]
    # record the total number of qrdr calls
    amr_out_dict['num QRDR'] = [num_qrdr_calls]
    # add a sulfamethoxazole column, only gets an R if there is both a dfr and a sul
    if sulfamethoxazole == 2:
        amr_out_dict['trimethoprim-sulfamethoxazole'] = ['R: ' + ';'.join(sulfamethoxazole_determinants)]
    else:
        amr_out_dict['trimethoprim-sulfamethoxazole'] = ['S']
    # make table, renaming final column to be the pST for the IncHI1 plasmid
    amr_out_df = pd.DataFrame(amr_out_dict)
    # rename column for plasmid ST
    amr_out_df.rename(columns={'NA_C19241A': 'IncHI1_ST6'}, inplace=True)
    return amr_out_df

def inspect_calls(full_lineage_data):
    # list of mykrobe genotypes to ignore
    ignore_genotype_levels = ['lineage1.2.0', 'lineage1.0', 'lineage1.2.3.0', 'lineage1.2.3.4.0', 'lineage1.2.3.4.0.3',
                              'lineage0', 'lineage0.0', 'lineage1.2.0.0', 'lineage1.2.3.0.0']
    # for all the genotypes that mykrobe is calling, inspect the values listed next to the actual call hierarhcy
    genotype_details = full_lineage_data['calls_summary']
    genotype_list = list(genotype_details.keys())
    # intialise empty values for best score and corresponding genotype
    best_score = 0
    best_genotype = None
    # node support is the number of 'good nodes' called by mykrobe
    # needs to be X/Y, where Y is the total number of levels at that genotype
    node_support = None
    # loop through each genotype
    for genotype in genotype_list:
        # the maximum score we can get is the total depth of the heirarchy
        max_score = genotype_details[genotype]['tree_depth']
        # the actual score is a sum of the values within each levels call
        actual_score = sum(list(genotype_details[genotype]['genotypes'].values()))
        # if actual score is equal to the max score, then this is the best genotype
        # BUT WE NEED TO DEAL WITH INSTANCES WHERE WE MIGHT HAVE A GREAT CALL 3.7.25 say AND A LESS GREAT CALL 3.6.1.1.2 say BUT BECACUSE THE HEIRARHCY IS BIGGER FOR 3.6.1.1.2 WE INADVERTANTLY CALL THAT
        if actual_score == max_score:
            best_score = actual_score
            best_genotype = genotype
            #node_support = genotype_details[best_genotype]['good_nodes']
        # if actual score is < max score, but is greater than the current best score,
        # then this is our best genotype for the moment
        elif actual_score < max_score and actual_score > best_score:
            best_score = actual_score
            best_genotype = genotype
            #node_support = genotype_details[best_genotype]['good_nodes']

    # get the best genotype in the right format
    best_genotyphi_genotype = best_genotype

    # determine any quality issues and split them amongst our various columns
    # if the node support for the best genotype is not '1' at all positions, we need to report that
    best_calls = genotype_details[best_genotype]['genotypes']
    # remove any calls that are prepended with "lineage", as these are fake levels in the hierarchy
    best_calls = dict([(x, y) for x, y in best_calls.items() if not x.startswith("lineage")])
    # get a list of all the calls for calculating confidence later
    best_calls_vals = list(best_calls.values())
    poorly_supported_markers = []
    # dict that keeps percentage supports, key=percent support; value=level
    lowest_within_genotype_percents = {}
    # list that gives supports for all markers in best call
    final_markers = []
    for level in best_calls.keys():
      # regardless of the call, get info
      call_details = full_lineage_data['calls'][best_genotype][level]
      # check that there is something there
      if call_details:
          # need to do this weird thing to grab the info without knowing the key name
          call_details = call_details[list(call_details.keys())[0]]
          ref = call_details['info']['coverage']['reference']['median_depth']
          alt = call_details['info']['coverage']['alternate']['median_depth']
          # calculate percent support for marker
          try:
              percent_support = alt / (alt + ref)
          except ZeroDivisionError:
              percent_support = 0
          # get the level name with the genotyphi call, ignoring any markers that don't exist
          if not level.startswith("lineage"):
              marker_string = level + ' (' + str(best_calls[level]) + '; ' + str(alt) + '/' + str(ref) + ')'
              final_markers.append(marker_string)
      # if the value is null, just report 0 (indicates that no SNV detected, either ref or alt?)
      else:
          lowest_within_genotype_percents[0] = level
          # get the level name with the genotyphi call, ignoring any markers that don't exist
          if not level.startswith("lineage"):
              marker_string = level + ' (0)'
              # its returned null so is therefore by definition poorly supported
              poorly_supported_markers.append(marker_string)
              # add to the final markers list too
              final_markers.append(marker_string)
              percent_support = 0
      # if call is 1 then that is fine, count to determine confidence later
      # if call is 0.5, then get info
      # if call is 0, there will be no info in the calls section, so just report 0s everywhere
      if best_calls[level] < 1:
          # then it must be a 0 or a 0.5
          # report the value (0/0.5), and also the depth compared to the reference
          # get the level name with the genotyphi call, ignoring any markers that don't exist
          if not level.startswith("lineage"):
              lowest_within_genotype_percents[percent_support] = level
              # add this to the list of poorly supported markers
              poorly_supported_markers.append(marker_string)

    # determining final confidence is based ONLY on the actual genotype, not incongruent genotype calls
    # strong = quality 1 for all calls
    if best_calls_vals.count(0) == 0 and best_calls_vals.count(0.5) == 0:
        confidence = 'strong'
        lowest_support_val = '-'
    # moderate = quality 1 for all calls bar 1 (which must be quality 0.5 with percent support > 0.5)
    elif best_calls_vals.count(0) == 0 and best_calls_vals.count(0.5) == 1:
        if min(lowest_within_genotype_percents.keys()) > 0.5:
            confidence = 'moderate'
        else:
            confidence = 'weak'
        lowest_support_val = round(min(lowest_within_genotype_percents.keys()), 3)
    # weak = more than one quality < 1, or the single 0.5 call is < 0.5% support
    else:
        confidence = 'weak'
        lowest_support_val = round(min(lowest_within_genotype_percents.keys()), 3)

    # make a list of all possible quality issues (incongruent markers, or not confident calls within the best geno)
    non_matching_markers = []
    non_matching_supports = []
    # we now want to report any additional markers that aren't congruent with our best genotype
    #ie if 3.6.1 is the best genotype, but we also have a 3.7.29 call, we need to report the 3.7 and 3.29 markers as incongruent
    if len(genotype_list) > 1:
        # remove the best genotype from the list
        genotype_list.remove(best_genotype)
        # loop through each incongruent phenotype
        for genotype in genotype_list:
            # extract the calls for that genotype
            other_calls = genotype_details[genotype]['genotypes']
            # for every call in our BEST calls, we're only interested
            # in calls that are incongruent
            for call in other_calls.keys():
                if call not in best_calls.keys():
                    call_info = full_lineage_data['calls'][genotype][call]
                    # check that there is something there (if the value is null, don't report it for incongruent calls)
                    if call_info:
                        # need to do this weird thing to grab the info without knowing the key name
                        call_info = call_info[list(call_info.keys())[0]]
                        ref_depth = call_info['info']['coverage']['reference']['median_depth']
                        alt_depth = call_info['info']['coverage']['alternate']['median_depth']
                        # only keep the call if the alternate has a depth of > 1
                        # this is because mykrobe fills in intermediate levels of the hierarchy with 0s
                        # if a lower level SNV marker is detected
                        if alt_depth >= 1:
                            percent_support = alt_depth / (alt_depth + ref_depth)
                            non_matching_supports.append(percent_support)
                            # get genotyphi call
                            if not call.startswith("lineage"):
                                marker_string = call + ' (' + str(other_calls[call]) + '; ' + str(alt_depth) + '/' + str(ref_depth) + ')'
                                non_matching_markers.append(marker_string)

    # get max value for non matching supports, if empty, return ''
    if len(non_matching_supports) > 0:
        max_non_matching = round(max(non_matching_supports), 3)
    else:
        max_non_matching = '-'
    # add '-' for those columns that are empty
    if len(poorly_supported_markers) == 0:
        poorly_supported_markers = ['-']
    if len(non_matching_markers) == 0:
        non_matching_markers = ['-']
    if len(final_markers) == 0:
        final_makers = ['-']

    return best_genotyphi_genotype, confidence, lowest_support_val, poorly_supported_markers, max_non_matching, non_matching_markers, final_markers

def extract_lineage_info(lineage_data, genome_name):

    # first, extract phylo spp information - if the spp is not Typhi then don't proceed
    spp_data = lineage_data['species']
    spp_call = list(spp_data.keys())[0]
    column_order = ['genome', 'species', 'spp_percent', 'genotype', 'confidence',
                    'lowest support for genotype marker', 'poorly supported markers',
                    'node support', 'max support for additional markers',
                    'additional markers']
    # if spp is unknown, then this is not typhi, exit this function
    if spp_call == "Unknown":
        out_dict = {'genome': [genome_name], 'species': ['not typhi'], 'spp_percent': [0], 'genotype': ['NA'],
                    'confidence': ['NA'], 'lowest support for genotype marker': ['NA'], 'poorly supported markers': ['NA'],
                    'node support': ['NA'], 'max support for additional markers': ['NA'], 'additional markers': ['NA']}
        out_df = pd.DataFrame(out_dict, columns=column_order)
        return out_df, spp_call
    else:
        # if it is typhi, then get the percentage
        spp_percentage = spp_data["Salmonella_Typhi"]["percent_coverage"]
        # if the percentage is <89, then exit this function as it's likely not typhi
        if spp_percentage < 85:
            out_dict = {'genome':[genome_name], 'species':['not typhi'], 'spp_percent': [spp_percentage],
                        'genotype': ['NA'], 'confidence': ['NA'], 'lowest support for genotype marker': ['NA'],
                        'poorly supported markers': ['NA'], 'node support': ['NA'], 'max support for additional markers': ['NA'],
                        'additional markers':['NA']}
            out_df = pd.DataFrame(out_dict, columns=column_order)
            return out_df, "Unknown"

    # if we are typhi, then continue
    #set up dictionary for final table output
    lineage_out_dict = {'genome': [genome_name]}

    # this try/except statement deals with instances where the genome is typhi but no markers are detected
    # in this instace return 'lineage0'
    try:
        genotype_calls = lineage_data['lineage']['lineage']
    except KeyError:
        out_dict = {'genome': [genome_name], 'species': ['typhi'], 'spp_percent': [spp_percentage],
                    'genotype': ['lineage0'], 'confidence': ['NA'], 'lowest support for genotype marker': ['NA'],
                    'poorly supported markers': ['NA'], 'max support for additional markers': ['NA'],
                    'additional markers': ['NA'], 'node support': ['NA']}
        out_df = pd.DataFrame(out_dict, columns=column_order)
        return out_df, spp_call
    # if there are no calls, populate with none
    if len(genotype_calls) == 0:
        lineage_out_dict['final_genotype'] = ['uncalled']
        lineage_out_dict['confidence'] = ['NA']
        lineage_out_dict['lowest support for genotype marker'] = ['NA']
        lineage_out_dict['poorly supported markers'] = ['NA']
        lineage_out_dict['max support for additional markers'] = ['NA']
        lineage_out_dict['additional markers'] = ['NA']
        lineage_out_dict['node support'] = ['NA']
    # if there is just one call, list that - but CHECK that all values in heirarchy are good (>0.5)
    # only report up to level in heirarchy where we have a good call
    if len(genotype_calls) == 1:
        # inspect calls for genotype
        best_genotype, confidence, lowest_support_val, poorly_supported_markers, non_matching_support, non_matching_markers, final_markers = inspect_calls(lineage_data['lineage'])
        # extract genotype from lineage thing, add that to the table
        #genotype = best_genotype.split('lineage')[1]
        lineage_out_dict['genotype'] = [best_genotype]
        # fill in all details on confidence/support
        lineage_out_dict['confidence'] = [confidence]
        lineage_out_dict['lowest support for genotype marker'] = [lowest_support_val]
        lineage_out_dict['poorly supported markers'] = ['; '.join(poorly_supported_markers)]
        lineage_out_dict['max support for additional markers'] = [non_matching_support]
        lineage_out_dict['additional markers'] = ['; '.join(non_matching_markers)]
        lineage_out_dict['node support'] = ['; '.join(final_markers)]
        #lineage_out_dict['all_genotype_calls'] = genotype_calls
    # if there is more than one call, we want to report the best, but also other calls
    elif len(genotype_calls) > 1:
        # get the info for each call
        # work out which lineage has the best call
        best_genotype, confidence, lowest_support_val, poorly_supported_markers, non_matching_support, non_matching_markers, final_markers = inspect_calls(lineage_data['lineage'])
        # now write out info
        #lineage_out_dict['genotype'] = [best_genotype.split('lineage')[1]]
        lineage_out_dict['genotype'] = [best_genotype]
        lineage_out_dict['confidence'] = [confidence]
        lineage_out_dict['lowest support for genotype marker'] = [lowest_support_val]
        lineage_out_dict['poorly supported markers'] = ['; '.join(poorly_supported_markers)]
        lineage_out_dict['max support for additional markers'] = [non_matching_support]
        lineage_out_dict['additional markers'] = ['; '.join(non_matching_markers)]
        lineage_out_dict['node support'] = ['; '.join(final_markers)]
    # add species info
    lineage_out_dict['species'] = ['typhi']
    lineage_out_dict['spp_percent'] = [spp_percentage]
    lineage_out_df = pd.DataFrame(lineage_out_dict, columns=column_order)

    return lineage_out_df, spp_call

def main():

    args = get_arguments()
    # list of tables for each result type
    results_tables = []

    # read in json files
    for json_file in args.jsons:
        with open(json_file) as f:
            myk_result = json.load(f)
        # get genome name (should first and ONLY key at start of json file)
        if len(list(myk_result.keys())) > 1:
            print("More than one result in mykrobe output file " + json_file + ", quitting")
            sys.exit()
        genome_name = list(myk_result.keys())[0]
        print(genome_name)
        # extract all the data for this genome
        genome_data = myk_result[genome_name]
        # extract the genotyping information
        lineage_data = genome_data["phylogenetics"]
        lineage_table, spp_call = extract_lineage_info(lineage_data, genome_name)
        # combine the two tables together
        # extract the amr information
        genome_amr_table = extract_amr_info(genome_data, genome_name, spp_call)
        # combine the two tables together
        all_info_table = pd.merge(lineage_table, genome_amr_table, on="genome", how="inner")
        results_tables.append(all_info_table)

    # concatenate, re-order columns, and write out
    final_results = pd.concat(results_tables, sort=True)
    column_order = ["genome", "species", "spp_percent", "genotype", "confidence",
                    "lowest support for genotype marker", "poorly supported markers",
                    "max support for additional markers", "additional markers", "node support", "ampicillin",
                    "azithromycin", "ceftriaxone", "ciprofloxacin", "chloramphenicol", "sulfonamides",
                    "trimethoprim", "trimethoprim-sulfamethoxazole", "tetracycline", 'IncFIAHI1', 'IncHI1A', 'IncHI1BR27', 'IncHI1_ST6', 'IncY', 'IncX3',
                    'IncHI2A', 'IncI1', 'IncL_M', 'IncFIB_pHCM2', 'IncFIB_K', 'IncN', 'z66', 'num QRDR',
                    'parC_S80I', 'parC_S80R', 'parC_E84G', 'parC_E84K', 'gyrA_S83F', 'gyrA_S83Y', 'gyrA_D87G', 'gyrA_D87N',
                    'gyrA_D87V', 'gyrA_D87Y', 'gyrB_S464F', 'gyrB_S464Y', 'acrB_R717L', 'acrB_R717Q',
                    'mphA', 'ermB', 'ereA', 'blaTEM-1D', 'blaCTX-M-15', 'AmpC1', 'blaOXA-7', 'blaOXA-134', 'blaSHV-12', 'qnrS1', 'qnrB1', 'qnrD1', 'catA1', 'cmlA1',
                    'sul1', 'sul2', 'dfrA1', 'dfrA5', 'dfrA7', 'dfrA14', 'dfrA15', 'dfrA17', 'dfrA18', 'tetA',
                    'tetB', 'tetC', 'tetD']
    final_results = final_results.reindex(columns=column_order)
    final_results.to_csv(args.prefix + "_predictResults.tsv", index=False, sep="\t")


if __name__ == '__main__':
    main()
