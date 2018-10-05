cd /media/Space2/home/andreas/miRNAschisto/sequences
awk '{print $1}' NCBI_11140_proteins_Shae.fasta | fa2tab | awk -F"\t" 'BEGIN{while(getline < "11140_lookup_IDs_NCBI.tsv"){look[">"$1]=$2}}{if($1 in look){print ">"look[$1]; print $2}else{print "ERROR"}}' > MS3_NCBI_11140_ol.fasta
mv MS3_NCBI_11140_ol.fasta ~/DrugTargetBaseLineShae/
cd ~/DrugTargetBaseLineShae/
mv MS3_NCBI_11140_ol.fasta reanalysis
cd reanalysis
grep ">" MS3_NCBI_11140_ol.fasta | sed 's/>//' > MS3_NCBI_11140_ol.ids
blastp -db ../MmBLASTs/Mus_musculus.GRCm38.pep.all.fa -query MS3_NCBI_11140_ol.fasta  -out blastpShMm.tsv -evalue 1e-5 -outfmt 6 -num_threads 12 2>&1 > /dev/null
blastp -db ../CeBLASTs/c_elegans.PRJNA13758.WS262.protein.fa -query MS3_NCBI_11140_ol.fasta  -out blastpShCe.tsv -evalue 1e-5 -outfmt 6 -num_threads 12 2>&1 > /dev/null
blastp -db ../DmBLASTs/dmel-all-translation-r6.19.fasta -query MS3_NCBI_11140_ol.fasta -out blastpShDm.tsv -evalue 1e-5 -outfmt 6 -num_threads 12 2>&1 > /dev/null

bestBLAST blastpShMm.tsv
awk -F"\t" 'BEGIN{while(getline < "../MmBLASTs/lethal_Mm_proteins.ids"){leth[$0]=$0}}{split($2,prot,/\./); if(prot[1] in leth){print $1"\t"$2}}' best_hits_blastpShMm.tsv > lethal_ShMm.tsv
bestBLAST blastpShCe.tsv
awk -F"\t" 'BEGIN{while(getline < "../CeBLASTs/lethal_proteins.ls"){leth[$0]=$0}}{if($2 in leth){print $1"\t"$2}}' best_hits_blastpShCe.tsv > lethal_ShCe.tsv
bestBLAST blastpShDm.tsv
awk -F"\t" 'BEGIN{while(getline < "../DmBLASTs/lethal_protein_ids.tsv"){leth[$0]=$0}}{if($2 in leth){print $1"\t"$2}}' best_hits_blastpShDm.tsv > lethal_ShDm.tsv
awk -F"\t" 'BEGIN{while(getline < "lethal_ShMm.tsv"){mm[$1]=$2}; while(getline < "lethal_ShCe.tsv"){ce[$1]=$2}; while(getline < "lethal_ShDm.tsv"){dm[$1]=$2}; print "#MS3_id\tcel_leth\tdmel_leth\tmus_leth"}{OFS="\t"; print $0,ce[$0],dm[$0],mm[$0]}' MS3_NCBI_11140_ol.ids > lethality_Ce_Dm_Mm_final.tsv
blastp -db /media/Space1/Blastdb/SwissProt2017May/uniprot_sprotdb -query MS3_NCBI_11140_ol.fasta  -out blastpShSwPro.tsv -evalue 1e-5 -outfmt 6 -num_threads 12 2>&1 > /dev/null
blastp -db /media/Space1/Blastdb/KEGG2017May/kegg_eukaryotesdb -query MS3_NCBI_11140_ol.fasta  -out blastpShKEGG.tsv -evalue 1e-5 -outfmt 6 -num_threads 12 2>&1 > /dev/null
blastp -db ../HsBLASTs/Homo_sapiens.GRCh38.pep.all.fa -query MS3_NCBI_11140_ol.fasta  -out blastpShHs.tsv -evalue 1e-5 -outfmt 6 -num_threads 12 2>&1 > /dev/null
cut -f1,2 best_hits_blastpShCe.tsv > best_hits_blastpShCe.ids
cut -f1,2 best_hits_blastpShMm.tsv > best_hits_blastpShMm.ids
cut -f1,2 best_hits_blastpShDm.tsv > best_hits_blastpShDm.ids
needle_new.sh best_hits_blastpShCe.ids ../data/nodots_Schistosoma_haematobium_v2.proteins.final.fa c_elegans.PRJNA13758.WS262.protein.fa 2>&1 > /dev/null
awk -F"\t" 'BEGIN{while(getline < "../data/gene_convert_table_update_15SEP2015.txt"){conv[$1]=$2}; while(getline < "../data/RNA-Seq/RSEM_female.genes.results"){female[$1]=$6}; while(getline < "../data/RNA-Seq/RSEM_male.genes.results"){male[$1]=$6}; print "#MS3_id\tmale_tpm\tfemale_tpm"}{if($0 in conv){if(conv[$0] in male){printf $0"\t"male[conv[$0]]}else{printf $0"\t0.00"}; if(conv[$0] in female){printf "\t"female[conv[$0]]}else{printf "\t0.00"}; print""}else{print $0"\t0.00\t0.00"}}' MS3_NCBI_11140_ol.ids > transcription_m_f_final.tsv
bestBLAST blastpShHs.tsv
cut -f1,2 best_hits_blastpShHs.tsv > best_hits_blastpShHs.ids
needle_new.sh best_hits_blastpShHs.ids MS3_NCBI_11140_ol.fasta ../HsBLASTs/Homo_sapiens.GRCh38.pep.all.fa 2>&1 > /dev/null
qtl75=$(cut -f4 output_best_hits_blastpShHs.ids | sort -n | awk -F"\t" '{arr[NR]=$0}END{for(i=1; i<=length(arr); i++){if(i/NR >= 0.75){print arr[i]; break}}}')
awk -F"\t" -v qtl=$qtl75 'OFS="\t"{$6=$5; if($4 <= qtl){$5="TRUE"}else{$5="FALSE"}; print $0}' output_best_hits_blastpShHs.ids > tmp
mv tmp output_best_hits_blastpShHs.ids
needle_new.sh best_hits_blastpShMm.ids MS3_NCBI_11140_ol.fasta ../MmBLASTs/Mus_musculus.GRCm38.pep.all.fa 2>&1 > /dev/null
needle_new.sh best_hits_blastpShDm.ids MS3_NCBI_11140_ol.fasta ../DmBLASTs/dmel-all-translation-r6.19.fasta 2>&1 > /dev/null
bestBLAST blastpShSwPro.tsv
cut -f1,2 best_hits_blastpShSwPro.tsv > best_hits_blastpShSwPro.ids

needle_new.sh best_hits_blastpShSwPro.ids MS3_NCBI_11140_ol.fasta MS3_NCBI_11140_ol.fasta /media/Space1/Blastdb/SwissProt2017May/uniprot_sprot.fasta 2>&1 > /dev/null
keggAnno -i blastpShKEGG.tsv -d keggAnno 2>&1 > /dev/null
grep -v "^##" ./keggAnno/blastpShKEGG.tsv.kegg.pathway.details | awk -F"\t" 'BEGIN{while(getline < "./keggAnno/blastpShKEGG.tsv.kegg.pathway.details"){pw[$2";"$4]++;}}{if(pw[$2";"$4]==1){print $1"\t"$2"\t"$4}}' | cut -f1 | sort | uniq > chokepoints.ls
ls ../flukesBLASTs/*.WBPS9.protein.fa | while read line; do blastp -db $line -query MS3_NCBI_11140_ol.fasta -out blastp$line.tsv -evalue 1e-5 -outfmt 6 -num_threads 12; done
ls blastp*.WBPS9.protein.fa.tsv | while read line; do bestBLAST $line; done
cat ../flukesBLASTs/*.WBPS9.protein.fa > all_flukes_protein.fa
ls best_hits*.WBPS9.protein.fa.tsv | while read line; do cut -f1,2 $line > "$line"_tmp; needle_new.sh "$line"_tmp MS3_NCBI_11140_ol.fasta all_flukes_protein.fa; rm "$line"_tmp; done
psiblast -db ../DrugBank/clean_drugbank.fasta -query MS3_NCBI_11140_ol.fasta -out psiblast_DruBa.tsv -evalue 1e-30 -outfmt "6 std stitle" -num_iterations 3 2>&1 > /dev/null
sed -i '/^$/d' psiblast_DruBa.tsv
sed -i '/Search has CONVERGED\!/d' psiblast_DruBa.tsv
bestBLAST psiblast_DruBa.tsv
awk -F"\t" 'BEGIN{while(getline < "../DrugBank/clean_drugbank.fasta"){if($0 ~ /^>/){prot[$1]=$2}}; while(getline < "DrugBank5.0.11.info"){db[$1]=$0}}{split(prot[">"$2],drugs,/;/); for(i in drugs){if(drugs[i] in db){print $1"\t"$2"\t"db[drugs[i]]}}}' best_hits_psiblast_DruBa.tsv > Sh_DB5_target_cmpd_assoc.tsv
psiblast -db ../ChEMBL/filtered_chembl23.fasta -query MS3_NCBI_11140_ol.fasta -out psiblast_chembl23.tsv -evalue 1e-30 -outfmt "6 std stitle" -num_iterations 3 2>&1 > /dev/null
sed -i '/^$/d' psiblast_chembl23.tsv
sed -i '/Search has CONVERGED\!/d' psiblast_chembl23.tsv
bestBLAST psiblast_chembl23.tsv
cut -f2 best_hits_psiblast_chembl23.tsv | sort | uniq > accession_targets.ids

#LAST STEPS
#on leonard!
./final_chemblSql23.py accession_targets.ids > accession_targets_chembl.tsv
sed -i 's/\t$//' accession_targets_chembl.tsv

#awk -F"\t" 'BEGIN{while(getline < "lookup_MS3_acc.tsv"){MS3[$1]=$2}}{for(i in MS3){if(MS3[i]==$1){print i"\t"$0}}}' 1133_accession_targets_chembl.tsv > final_1133_accession_targets_chembl.tsv
#awk -F"\t" '{if($5 == 0 || $5 == "None"){print $0}}' final_1133_accession_targets_chembl.tsv > ro5_filtered.tsv
#awk -F"\t" '{if($17 == "B" && ($18 == "Activity" || $18 == "Kd" || $18 == "Ki" || $18 == "Inhibition" || $18 == "Potency" || $18 == "IC50")){print $0}}' ro5_filtered.tsv > activity_type_filtered.tsv
awk -F"\t" '{if($17 == "B" && ($18 == "Activity" || $18 == "Kd" || $18 == "Ki" || $18 == "Inhibition" || $18 == "Potency" || $18 == "IC50")){print $0}}' accession_targets_chembl.tsv > noro5_activity_type_filtered.tsv
awk -F"\t" '{if((($18 == "Activity" || $18 == "Inhibition") && ($20 == "%")) || (($18 == "Kd" || $18 == "Potency" || $18 == "Ki" || $18 == "IC50") && ($20 == "nM"))){print $0}}' noro5_activity_type_filtered.tsv | awk -F"\t" '{if($19 == "None" || $21 == "None"){}else{print $0}}' > noro5_unit_relation_filtered.tsv
#awk -F"\t" '{if((($18 == "Activity" || $18 == "Inhibition") && ($20 == "%")) || (($18 == "Kd" || $18 == "Potency" || $18 == "Ki" || $18 == "IC50") && ($20 == "nM"))){print $0}}' activity_type_filtered.tsv | awk -F"\t" '{if($19 == "None" || $21 == "None"){}else{print $0}}' > unit_relation_filtered.tsv
#awk -F"\t" '{if(($20 == "nM" && $21 <= 10000) || ($20 == "%" && $21 >= 70)){print $0}}' unit_relation_filtered.tsv > value_filtered.tsv
awk -F"\t" '{if(($20 == "nM" && $21 <= 10000) || ($20 == "%" && $21 >= 70)){print $0}}' noro5_unit_relation_filtered.tsv > noro5_value_filtered.tsv



#compile final table
awk -F"\t" 'BEGIN{while(getline < "MS3_NCBI_11140_ol.ids"){look[$0]=$0}; while(getline < "output_best_hits_schistosoma_mansoni.PRJEA36577.WBPS9.protein.fa.tsv_blastp_tmp"){sman[$1]=$2"\t"$3"\t"$4"\t"$5}; for(i in look){if(sman[i]==""){sman[i]="\t\t\t"}}; while(getline < "output_best_hits_schistosoma_japonicum.PRJEA34885.WBPS9.protein.fa.tsv_blastp_tmp"){sjap[$1]=$2"\t"$3"\t"$4"\t"$5}; for(i in look){if(sjap[i]==""){sjap[i]="\t\t\t"}}; while(getline < "output_best_hits_opisthorchis_viverrini.PRJNA222628.WBPS9.protein.fa.tsv_blastp_tmp"){opvi[$1]=$2"\t"$3"\t"$4"\t"$5}; for(i in look){if(opvi[i]==""){opvi[i]="\t\t\t"}}; while(getline < "output_best_hits_clonorchis_sinensis.PRJDA72781.WBPS9.protein.fa.tsv_blastp_tmp"){closi[$1]=$2"\t"$3"\t"$4"\t"$5}; for(i in look){if(closi[i]==""){closi[i]="\t\t\t"}}; while(getline < "output_best_hits_fasciola_hepatica.PRJNA179522.WBPS9.protein.fa.tsv_blastp_tmp"){fahe[$1]=$2"\t"$3"\t"$4"\t"$5}; for(i in look){if(fahe[i]==""){fahe[i]="\t\t\t"}}; print "MS3_id\tsman_id\tsman_ident\tsman_simil\tsman_cov\tsjap_id\tsjap_ident\tsjap_simil\tsjap_cov\topvi_id\topvi_ident\topvi_simil\topvi_cov\tclosi_id\tclosi_ident\tclosi_simil\tclosi_cov\tfahe_id\tfahe_ident\tfahe_simil\tfahe_cov"}{OFS="\t"; print $0,sman[$0],sjap[$0],opvi[$0],closi[$0],fahe[$0]}' MS3_NCBI_11140_ol.ids > fluke_BLASTs_final.tsv

awk -F"\t" 'BEGIN{while(getline < "MS3_NCBI_11140_ol.ids"){look[$0]=$0}; while(getline < "output_best_hits_blastpShCe.ids"){cel[$1]=$2"\t"$3"\t"$4"\t"$5}; for(i in look){if(cel[i]==""){cel[i]="\t\t\t"}}; while(getline < "output_best_hits_blastpShDm.ids"){dmel[$1]=$2"\t"$3"\t"$4"\t"$5}; for(i in look){if(dmel[i]==""){dmel[i]="\t\t\t"}}; while(getline < "output_best_hits_blastpShMm.ids"){mus[$1]=$2"\t"$3"\t"$4"\t"$5}; for(i in look){if(mus[i]==""){mus[i]="\t\t\t"}}; while(getline < "output_best_hits_blastpShHs.ids"){hsa[$1]=$2"\t"$3"\t"$4"\t"$5"\t"$6}; for(i in look){if(hsa[i]==""){hsa[i]="\t\t\t\t"}}; print "MS3_id\tcel_id\tcel_ident\tcel_simil\tcel_cov\tdmel_id\tdmel_ident\tdmel_simil\tdmel_cov\tmus_id\tmus_ident\tmus_simil\tmus_cov\thsa_id\thsa_ident\thsa_simil\thsa_lower_75_qtl\thsa_cov"}{OFS="\t"; print $0,cel[$0],dmel[$0],mus[$0],hsa[$0]}' MS3_NCBI_11140_ol.ids > other_species_BLASTs.tsv

awk -F"\t" 'BEGIN{while(getline < "chokepoints.ls"){choke[$0]="TRUE"}; print "MS3_id\tchokepoint"}{if($0 in choke){print $0"\t"choke[$0]}else{print $0"\tFALSE"}}' MS3_NCBI_11140_ol.ids > chokepoints_final.tsv

awk -F"\t" 'BEGIN{print "MS3_id\ttarget_id\tcmpd_id\tcmpd_name\tstatus\tsmiles\tro5_pass\tmddr-like"}{print $0}' Sh_DB5_target_cmpd_assoc.tsv > DrugBank_Sh_final.tsv

awk -F"\t" 'BEGIN{while(getline < "fluke_BLASTs_final.tsv"){fluke[$1]=$0}; while(getline < "transcription_m_f_final.tsv"){trans[$1]=$0}; while(getline < "other_species_BLASTs.tsv"){other[$1]=$0}; while(getline < "lethality_Ce_Dm_Mm_final.tsv"){let[$1]=$0}; while(getline < "chokepoints_final.tsv"){choke[$1]=$0}; OFS="\t"; for(i in fluke){print fluke[i],trans[i],other[i],let[i],choke[i]}}' > tmp
cut -f1-21,23-24,26-42,44-46,48 tmp > tmp2

sort -r tmp2 > merged_Sh_ranking_final.tsv


sed 's/|/\t/g' output_best_hits_blastpShSwPro.ids | cut -f1,3- | awk -F"\t" 'BEGIN{print "MS3_id\tuniprot_id\tgene_name_species\tsp_ident\tsp_sim\tsp_cov"}{print $0}' > tmp

awk -F"\t" 'BEGIN{while(getline < "tmp"){look[$1]=$0}}{if($1 in look){print look[$1]"\t"$0}else{print $1"\t\t\t\t\t\t"$0}}' merged_Sh_ranking_final.tsv > tmp2
cut -f1-6,8- tmp2 | sort -r > merged_Sh_ranking_final.tsv
