#!/usr/bin/env python3
import argparse
import os
import sys
import re


def argument_check(args):
    """
    Process and check user inputs
    """
    if os.path.isfile(args.lineage_file[0]) is False:
        sys.stderr.write("GTDB lineage file ({}) not found\n".format(args.lineage_file))
        sys.exit(1)
    if os.path.isdir(args.ncbi_taxdmp_dir[0]) is False:
        sys.stderr.write("NCBI taxdmp_dir ({}) was not found\n".format(args.taxdmp_dir))
        sys.exit(1)
    if os.path.isfile("{}/nodes.dmp".format(args.ncbi_taxdmp_dir[0])) is False:
        sys.stderr.write("NCBI taxdmp node file ({}/nodes.dmp) was not found\n".format(args.taxdmp_dir))
        sys.exit(1)
    if os.path.isfile("{}/names.dmp".format(args.ncbi_taxdmp_dir[0])) is False:
        sys.stderr.write("NCBI taxdmp name file ({}/names.dmp) was not found\n".format(args.taxdmp_dir))
        sys.exit(1)
    if os.path.dirname(args.outfile[0]) != "" and os.path.isdir(os.path.dirname(args.outfile[0])) is False:
        os.mkdir(os.path.dirname(args.outfile[0]))

def store_ncbi_txdmp(txdmp_dir):
    txdmp_dict = {"nodes": {}, "names_scientific": {}, "names_all": {}}
    with open("{}/nodes.dmp".format(txdmp_dir), 'r') as data:
        for line in data.readlines():
            line = re.sub("\t\|$", "", line.rstrip())
            split_line = re.split(r'\s\|\s', line.strip())
            txdmp_dict["nodes"][split_line[0]] = split_line[1:3]
            txdmp_dict["nodes"][split_line[0]][1] = txdmp_dict["nodes"][split_line[0]][1].replace("superkingdom", "domain")
    with open("{}/names.dmp".format(txdmp_dir), 'r') as data:
        for line in data.readlines():
            line = re.sub("\t\|$", "", line.rstrip())
            split_line = re.split(r'\s\|\s', line.strip())
            if split_line[3] in ["scientific name", "synonym", "equivalent name"]:
                txdmp_dict["names_all"][split_line[1]] = split_line[0]
                if split_line[3] == "scientific name":
                    txdmp_dict["names_scientific"][split_line[0]] = split_line[1]
    return txdmp_dict


def store_alternate_taxon_info(alternate_names_file):
    altnames_dict = {}
    with open(alternate_names_file) as data:
        for line in data.readlines():
            split_line = re.split(',', line.strip())
            altnames_dict[split_line[0]] = split_line[1:]
    return altnames_dict

def gtdb2ncbi(gtdb_lineage, txdmp_dict, altname_dict):
    gtax_list = gtdb_lineage.split(";")
    gtax_list.reverse()
    for gtax in gtax_list:
        gtax_split = re.split("__", gtax)
        ncbi_taxon = gtax_split[1]
        alt_name = False
        if gtax_split[1] in altname_dict.keys():
            ncbi_taxon = altname_dict[gtax_split[1]][0]
            ncbi_taxid = altname_dict[gtax_split[1]][1]
            alt_name = True
        else:
            if gtax_split[0] == "s" \
                    and gtax_split[1] not in txdmp_dict["names_all"].keys() \
                    and "Candidatus" + gtax_split[1] not in txdmp_dict["names_all"].keys():
                gtax_split[1] = gtax_split[1].split(" ")[0]
                gtax_split[0] = "g"
            if gtax_split[1] not in txdmp_dict["names_all"].keys() \
                    and "Candidatus" + gtax_split[1] in txdmp_dict["names_all"].keys():
                ncbi_taxon = "Candidatus" + gtax_split[1]
                gtax_split[1] = "Candidatus" + gtax_split[1]
            if gtax_split[0] not in ["s", "g"]:
                if gtax_list[-1] == "d__Archaea":
                    ncbi_taxon = gtax_split[1] + " archaeon"
                else:
                    ncbi_taxon = gtax_split[1] + " bacterium"
            elif gtax_split[0] == "g":
                ncbi_taxon = gtax_split[1] + " sp."
                if "uncultured " + gtax_split[1] in txdmp_dict["names_all"].keys():
                    ncbi_taxon = "uncultured " + gtax_split[1]
            if ncbi_taxon in txdmp_dict["names_all"].keys():
                ncbi_taxid = txdmp_dict["names_all"][ncbi_taxon]
                ncbi_taxon = txdmp_dict["names_scientific"][ncbi_taxid]
            else:
                ncbi_taxid = "NA"
        if ncbi_taxid != "NA":
            tmp_taxon = ncbi_taxon.replace(" bacterium", "").\
                replace(" archaeon", "").\
                replace(" cyanobacterium", "").\
                replace(" sp.", "").\
                replace("uncultured ", "")
            tmp_taxon = re.sub(' \(.*\)', '', tmp_taxon)
            if tmp_taxon not in txdmp_dict["names_all"].keys():
                tmp_taxon = tmp_taxon.title()
            tmp_taxid = txdmp_dict["names_all"][tmp_taxon]
            node = txdmp_dict["nodes"][tmp_taxid]
            ncbi_lineage = ""
            if node[1][0] != "s":
                ncbi_lineage = node[1][0] + "__" + tmp_taxon
            while node[1] != "domain":
                tmp_taxid = node[0]
                tmp_taxon = txdmp_dict["names_scientific"][tmp_taxid]
                node = txdmp_dict["nodes"][tmp_taxid]
                if node[1] in ["domain","kingdom","phylum","class","order","family","genus","species"]:
                    ncbi_lineage = node[1][0] + "__" + tmp_taxon + ";" + ncbi_lineage
            if ncbi_taxon == "bacterium":
                ncbi_taxon = "Bacteria bacterium"
            return ncbi_lineage + "," + ncbi_taxon + "," + ncbi_taxid + "\n"


def main(args):
    argument_check(args)
    txdmp_dict = store_ncbi_txdmp(args.ncbi_taxdmp_dir[0])
    altname_dict = store_alternate_taxon_info(args.alternate_names[0])
    unique_output_dict = {}
    output = open(args.outfile[0], "w")
    output.writelines("ncbi_classification,ncbi_taxon,taxon_id\n")
    with open(args.lineage_file[0], 'r') as data:
        for line in data.readlines():
            if line.strip() in unique_output_dict.keys():
                output.writelines(unique_output_dict[line.strip()])
            else:
                out = gtdb2ncbi(line.strip(),txdmp_dict, altname_dict)
                output.writelines(line.strip() + "," + out)
                unique_output_dict[line.strip()] = line.strip() + "," + out
    output.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Takes GTDB lineages and returns closest approximate NCBI lineages')
    parser.add_argument('lineage_file', metavar='lineage_file', type=str, nargs=1,
                        help='File containing GTDB lineages (one per line)')
    parser.add_argument('outfile', metavar='outfile', type=str, nargs=1,
                        help="Output CSV in format \"gtdb_lineage,ncbi_lineage,ncbi_taxon,ncbi_taxonID\"")
    parser.add_argument('-N', '--ncbi_taxdmp_dir', type=str, required=True, nargs=1,
                        help='NCBI taxdmp directory. Most recent data can be downloaded from '
                             '\'https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz\'')
    parser.add_argument('-A', '--alternate_names', type=str, nargs=1,
                        help='CSV of equivalent classifications not found in NCBI taxdmp files in the format '
                             '\"GTDB_taxon,NCBI_taxon,NCBI_taxonID\"')
    args = parser.parse_args()
    main(args)
