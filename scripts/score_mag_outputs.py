#!/usr/bin/env python3
import argparse
import os
import sys
import pandas as pd
import numpy as np

########################
# Quality Class Weighting
A = 0.5 # MEDIUM
B = 0.2 # Contamination
##########################

def sort_quals(qual):
    if qual == "HIGH":
        return(0)
    elif qual == "MEDIUM":
        return(1)
    else:
        return(2)

def format_data(bin_stats):
    # Import Data
    raw = pd.read_csv(bin_stats, sep=",")
    # Split name to get 'host', 'assembler', and 'part' columns
    split_name = raw.name.str.split(".", expand=True)
    split_name.columns = ["host", "assembler", "part"]
    # Put split names and data together
    bin_data = pd.concat([split_name, raw], axis=1)
    # Remove bins with greater than 10% contamination
    bin_data = bin_data[bin_data.Contamination <= 10.0]
    # Rename magscot and dastool outputs so they appear in binning program name
    for refiner in bin_data.refining_program.unique():
        bin_data.loc[bin_data['refining_program'] == refiner, 'binning_program'] = bin_data.loc[bin_data['refining_program'] == refiner, 'binning_program'].replace('.*.raw', "{}_raw".format(refiner), regex=True)
        bin_data.loc[bin_data['refining_program'] == refiner, 'binning_program'] = bin_data.loc[bin_data['refining_program'] == refiner, 'binning_program'].replace('.*.filtered', "{}_filtered".format(refiner), regex=True)
    # Remove any anomalous data
    bin_data = bin_data[(bin_data.contigs > 0)]
    bin_data = bin_data[~bin_data['part'].isin(['all_primary', 'snpt', 'alternate'])]
    bin_data.drop(['name'], axis=1, inplace=True)
    return(bin_data)

def score_bins(bin_data):
    bin_data['cch_score'] = (bin_data.Completeness - bin_data.Contamination - (bin_data['Strain heterogeneity']/100 * bin_data.Contamination))
    bin_data['contiguity_score'] = (bin_data.cch_score / bin_data.contigs)
    bin_data['mg_score'] = 0
    bin_data.loc[bin_data['unique_trnas'] >= 18, 'mg_score'] += 25
    bin_data.loc[bin_data['rrna_23s'] == 'Y', 'mg_score'] += 25
    bin_data.loc[bin_data['rrna_16s'] == 'Y', 'mg_score'] += 25
    bin_data.loc[bin_data['rrna_5s'] == 'Y', 'mg_score'] += 25
    bin_data.loc[bin_data['unique_trnas'] < 18, 'mg_score'] += np.round(25 * (bin_data.loc[bin_data['unique_trnas'] < 18, 'unique_trnas']/ 18))
    bin_data['score'] = (bin_data['cch_score'] + bin_data['contiguity_score'] + bin_data['mg_score'])/300
    bin_data.loc[bin_data['quality'] == 'MEDIUM', 'score'] *= A
    bin_data.loc[bin_data['quality'] == 'LOW', 'score'] *= B
    return(bin_data)
    
def create_count_table(bin_data):
    counts = bin_data[["host", "assembler", "part", "binning_program", "quality"]].value_counts().reset_index()
    scores = bin_data.groupby(["host", "assembler", "part", "binning_program",]).score.sum().reset_index()
    # Pivot the counts tables
    qual_values = bin_data.quality.unique().tolist()
    qual_values.sort()
    pivot_counts = counts.pivot(index=['host','assembler', 'part', 'binning_program'], columns='quality').reset_index().fillna(0)
    pivot_counts.columns = ["host", "assembler", "part", "binning_program"] + qual_values
    qual_values.sort(key=sort_quals)
    pivot_counts_reorder = pivot_counts[["host", "assembler", "part", "binning_program",] + qual_values]
    # Add columns of DREP assemblies and SSU counts to table
    ssu_counts = bin_data.groupby(["host", "assembler", "part", "binning_program"]).ssu_count.sum().reset_index()
    drep_counts = bin_data[bin_data.drep == "PASSED"][["host", "assembler", "part", "binning_program", "quality"]].value_counts().reset_index()
    drep_pivot_counts = drep_counts.pivot(index=['host','assembler', 'part', 'binning_program'], columns='quality').reset_index().fillna(0)
    drep_qual_values = bin_data[bin_data.drep == "PASSED"].quality.unique().tolist()
    drep_qual_values.sort()
    drep_pivot_counts.columns = ["host", "assembler", "part", "binning_program"] + drep_qual_values
    drep_qual_values.sort(key=sort_quals)
    drep_pivot_counts = drep_pivot_counts[["host", "assembler", "part", "binning_program"] + drep_qual_values]
    drep_pivot_counts.columns = ["host", "assembler", "part", "binning_program"] + [x + "_DREP" for x in drep_pivot_counts.columns[4:].tolist()]
    pivot_counts = pivot_counts_reorder.merge(drep_pivot_counts, how='outer')
    pivot_counts = pivot_counts.merge(ssu_counts)
    pivot_counts = pivot_counts.convert_dtypes()
    pivot_counts = pivot_counts.merge(scores)
    pivot_counts['mags'] = pivot_counts[qual_values].sum(axis=1)
    pivot_counts.sort_values(by=["score"], inplace=True, ascending=False)
    return(pivot_counts)
    
def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bin_stats", type=str, help="Path to bin data csv")
    parser.add_argument("-p", "--output_prefix", type=str, help="Prefix for output data", default="./map_out")
    args = parser.parse_args()
    if os.path.isfile(args.bin_stats) is False:
        sys.stderr.write("File ({}) not found\n".format(args.bin_stats))
        sys.exit(1)
    # Manipulate data
    bin_data = format_data(args.bin_stats)
    bin_data_scored = score_bins(bin_data)
    count_table = create_count_table(bin_data_scored)
    # Write output
    bin_data_scored = bin_data_scored.convert_dtypes()
    bin_data_scored.to_csv(args.output_prefix + ".scored_bindata.csv", index=False, na_rep="NA")
    count_table.to_csv(args.output_prefix + ".counts.csv", index=False, na_rep="NA")


main()
    
    
    
    
    