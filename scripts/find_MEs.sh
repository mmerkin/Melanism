datapath=BUSCO_GCA_905404145.2_ilBisBetu1.2_genomic.fna/run_lepidoptera_odb10/busco_sequences/single_copy_busco_sequences
tablepath=path/to/lookup/table
output_file=Biston_merian_elements.tsv


echo -e "busco_locus\tME\tchr\tstart\tstop" > "$output_file"

for file in "$datapath"/*.gff; do 
filetag=${file##*/}
filebase=${filetag%.gff}

match=$(grep -w "$filebase" busco/CW_merian_element_table.tsv)
[[ -z "$match" ]] && continue  
merian="${match#*$'\t'}"

mRNA=$(grep "mRNA" $datapath/$filetag)

IFS=$'\t' read -r chr _ _ start stop _ _ <<< "$mRNA"
echo -e "$filebase\t$merian\t$chr\t$start\t$stop" >> "$output_file"

done


#To do: take output, sort by second column, count number of uniq values for third column (can do in awk?), investigate suitable thresholds
