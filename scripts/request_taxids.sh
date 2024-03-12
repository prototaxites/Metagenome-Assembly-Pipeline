#!/usr/bin/env bash

# request_taxids.sh
# taxon_request_file format -> "proposed_name\tname_type\thost\tproject_id\tdescription"
# Name type options:
#	- "Environmental Name"
#	- "Sythetic Name"
#	- "Novel Species"
#	- "Unidentified Species"
#	- "Published Species"

ena_email="ena-asg@ebi.ac.uk"
sender_email="ng13@sanger.ac.uk"
taxon_request_file=$1
host_species="$2"

echo "Please review the attached taxonomy ID requests." \
	| mailx -s "Taxon ID requests for cobionts from the ASG species '$host_species'" \
		-A $taxon_request_file $ena_email $sender_email
