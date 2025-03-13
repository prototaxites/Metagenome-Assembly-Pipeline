#!/usr/bin/env bash

# metagenome_curation_request.sh

grit_email="grit-jira@sanger.ac.uk"
# grit_email="ng13@sanger.ac.uk"
yaml=$1

stats_line=`grep -n $'^stats:' $yaml | cut -d':' -f1`
name=`basename $yaml .yaml`
tail -n +$stats_line $yaml \
	| mailx -s "$name is ready for curation" \
		-A $yaml $grit_email
