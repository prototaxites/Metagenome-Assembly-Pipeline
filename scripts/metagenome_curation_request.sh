#!/usr/bin/env bash

# metagenome_curation_request.sh

grit_email="grit-jira@sanger.ac.uk"
# grit_email="ng13@sanger.ac.uk"
yaml=$1

name=`basename $yaml .yaml`
echo "" \
	| mailx -s "$name is ready for curation" \
		-A $yaml $grit_email
