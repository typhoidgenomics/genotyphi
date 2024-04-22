#!/usr/bin/env bash
set -ve
# make the background vcf file
python make_probes.make_vcf.py
# remove mongo-db folder
rm -rf $PWD/mongo-db

# initalise settings
ref_fa=AL513382.fasta
ref_gb=AL513382.gbk
db_name=typhi
db=$PWD/mongo-db/

# set up mongod
mkdir $db

mongod --quiet --dbpath $db &
sleep 5
mykrobe variants add -f --db_name $db_name make_probes.backgrounds.vcf $ref_fa -m CORTEX

# re-create Panel folder
rm -rf Panel
mkdir Panel

# make species probes
samtools faidx $ref_fa AL513382.1:2890400-2892457 | awk 'NR==1{$0 = ">invA?name=Salmonella_enterica&length=2057&panel_type=complex"} 1' | gzip -9 -c > Panel/typhi.probe.invA.20210322.fa.gz

# make MLST probes
./mlst_get_probes.py typhiSTs.txt pubmlst_salmonella Salmonella_Typhi tmp.mlst.probes
gzip -9 -c tmp.mlst.probes.fasta > Panel/typhi.mlst.20210322.fa.gz
rm -r tmp.mlst.probes.fasta tmp.mlst.probes.json 

# make qrdr and acrB probes
mykrobe variants make-probes --db_name $db_name -k21 -t qrdr.acrB.aa.in.txt -g $ref_gb $ref_fa > tmp.probes.aa.fa
# make plasmid pST probe
mykrobe variants make-probes --db_name $db_name -k21 -t pST_positions.txt --no-backgrounds AM412236.fasta > tmp.probes.pST.fa
#swap ref and alt in the pST probe
echo ">ref-C19241A?var_name=C19241A&num_alts=1&ref=AM412236&enum=0&gene=NA&mut=C19241A" > tmp.probes.reorder.pST.fa; sed -n 4p tmp.probes.pST.fa >> tmp.probes.reorder.pST.fa; echo ">alt-C19241A?var_name=C19241A&enum=0&gene=NA&mut=C19241A" >> tmp.probes.reorder.pST.fa; sed -n 2p tmp.probes.pST.fa >> tmp.probes.reorder.pST.fa
# make plasmid and AMR gene sequence probes
python make_probes.make_amr_pres_abs.py
# concatenate all probes together
cat tmp.probes.aa.fa typhi.amr.plas.presAbs.20240407.fa tmp.probes.reorder.pST.fa | gzip -9 -c > Panel/typhi.amr.plas.20240407.probes.fa.gz
# make json for amr and plasmid reps
python make_probes.make_res_to_drug_json.py > Panel/typhi.amr.plas.20240407.json

# make lineage probes
mykrobe variants make-probes --db_name $db_name -k21 -t typhi_lineage_positions.txt --lineage Panel/typhi.lineage.20221109.json $ref_fa | gzip -9 -c > Panel/typhi.lineage.20221109.probes.fa.gz

# make tarball
cp Panel/*.gz mykrobe.panel.typhi.20240407
cp Panel/*.json mykrobe.panel.typhi.20240407

tar -czvf mykrobe.panel.typhi.20240407.tar.gz mykrobe.panel.typhi.20240407

