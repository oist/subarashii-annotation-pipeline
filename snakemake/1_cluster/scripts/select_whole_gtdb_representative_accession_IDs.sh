#!/bin/bash

RELEASE=$1
OUTFILE=$2

wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release${RELEASE}/${RELEASE}.0/bac120_r${RELEASE}.tree
wget https://data.ace.uq.edu.au/public/gtdb/data/releases/release${RELEASE}/${RELEASE}.0/ar53_r${RELEASE}.tree

sed -E "s/'[0-9.]+\:[A-Za-z_0-9 ;-]+'//g" bac120_r${RELEASE}.tree > bac120_r${RELEASE}_fixed.tree
sed -E "s/'[0-9.]+\:[A-Za-z_0-9 ;-]+'//g" ar53_r${RELEASE}.tree > ar53_r${RELEASE}_fixed.tree

$(dirname "${BASH_SOURCE[0]}")/get_tree_leaves.py bac120_r${RELEASE}_fixed.tree > $OUTFILE
$(dirname "${BASH_SOURCE[0]}")/get_tree_leaves.py ar53_r${RELEASE}_fixed.tree >> $OUTFILE

rm bac120_r${RELEASE}_fixed.tree
rm bac120_r${RELEASE}.tree
rm ar53_r${RELEASE}_fixed.tree
rm ar53_r${RELEASE}.tree
