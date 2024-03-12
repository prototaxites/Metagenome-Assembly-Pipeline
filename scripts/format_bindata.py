#!/usr/bin/env python3

import pandas as pd
import argparse
import os

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("bin_stats", type=str, help="Path to bin data csv")
    parser.add_argument("-p", "--output_prefix", type=str, help="Prefix for output data", default="./out")
    args = parser.parse_args()
    if os.path.isfile(args.bin_stats) is False:
        sys.stderr.write("File ({}) not found\n".format(args.bin_stats))
        sys.exit(1)
    df = pd.read_csv(args.bin_stats, sep=",")
    df["length"] = df["size"].apply(lambda x: x / 10**6)
    outdf = df[["host","primary_tolid","primary_biosample","tolid","ncbi_taxon","taxon_id","bin_type","length","contigs","circular","Completeness","Contamination","mean_coverage","ssu_count","total_trnas","unique_trnas","rrna_23s","rrna_16s","rrna_5s","biosample"]]
    outdf = outdf.convert_dtypes()
    outdf.to_csv(args.output_prefix + ".csv", index=False, na_rep="NA")

main()
