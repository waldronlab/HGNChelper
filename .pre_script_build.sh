#!/bin/bash

GH_REPO="@github.com/waldronlab/HGNChelper.git"
FULL_REPO="https://$GH_TOKEN$GH_REPO"

pwd
cd inst/
Rscript --vanilla hgncLookup.R
cd ..

git init
git config user.name "lwaldron-travis"
git config user.email "travis"

git branch newmap
git checkout newmap
git stage data/hgnc.table.rda
git commit -m "updated HGNC map"
git push --force --quiet $FULL_REPO master:newmap

