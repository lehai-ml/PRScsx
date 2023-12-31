#!/bin/bash

#output_dir=test_data
#output_file=test
#ref_panels=ref_panels
#target_file=$src/$output_dir/test1
#summary_stats=$src/$output_file/EAS_sumstats_new.txt,$src/$output_file/EUR_sumstats_new.txt
#n_gwas=100000,200000
#populations=( EAS EUR )
#
while getopts "hb:c:l:t:o:n:p:d:" arg; do
    case $arg in
	b) IFS=, read -r -a base_files <<< "${OPTARG}";;
	c) IFS=/ read -r -a columns <<< "${OPTARG}" ;; #split the string by comma
	n) n_gwas=${OPTARG};;
	p) IFS=, read -r -a populations <<< "${OPTARG}" ;;
	l) ld_folder=${OPTARG};;
	t) target_file=${OPTARG};;
	d) output_dir=${OPTARG};;
	o) output_file=${OPTARG};;
	h) 
	    echo "-h: help menu"
	    echo "-b:path to summary statistics base files (must be more than 1). Make sure the file is not compressed .gz : Your summary statistics must be of the following format to work: SNP A1 A2 Beta(OR) P. If it is not in that format, you can provide the argument c "
	    echo "-c: this argument is used to call select the columns SNP A1 A2 Beta(OR) P from the base file provided in argument -b. If providing multiple files in base files, then columns for each corresponding file are sepearted by "/". eg. 1,2,3,3/1,2,3,4 For ASD(Grove2019) that is 2,4,5,7,9. SCZ(2022) is 2,4,5,9,11. ASD(Spark) is 2,5,6,7,9"
	    echo "NOTE: PRScs expect the SNP in 1000G and the other two files to match. So you may need to map chr:bp to rsid in all of your files. However, I have changed the PRScs.py to match by chr:bp instead of rsid."
	    echo "-p: populations of each respective base files. Separate by comma"
	    echo "-n: number of cases in the summary statistics. ASD(Grove) - 46350. SCZ (175799)"
	    echo "-l: path to the ld folder. This is downloaded from the original PRScs repository"
	    echo "-t: path to target b-file" 
	    echo "-d: path to output directory"
	    echo "-o: prefix of the output file"
	    exit 0;;
	?) exit 1;;
    esac
done
[[ -z "${base_files}" ]] && echo "ERROR: -b base files are missing" && exit 1
[[ -z "${target_file}" ]] && echo "ERROR: -t target file is missing" && exit 1
[[ -z "${ld_folder}" ]] && echo "ERROR: -l ld folder is missing" && exit 1
[[ -z "${populations}" ]] && echo "ERROR: -p populations names are missing " && exit 1
[[ -z "${output_dir}" ]] && echo "ERROR: -d output directory is missing " && exit 1
[[ -z "${output_file}" ]] && echo "ERROR: -o output file is missing " && exit 1
#
#the PRScs.py requires the SNPs in the reference (1000G) and the other two files to match. 
#So need to map chr:bp to rsid in all of your files.
#your summary statistics header must be in the following format. SNP A1 A2 BETA(or OR) P

# convenience function to output a new file tab separated table, by defining the column numbers.
# usage output_columns source_table output_table [columns separated by space]

function return_columns() {
    file=$1
    output=$2
    shift 2
    to_parse=($@)
    to_print=()
    for column in ${to_parse[@]}; do
        if [[ $column == ${to_parse[$(( ${#to_parse[*]} - 1 ))]} ]]; then
            to_print+=(\$${column})
        else
            to_print+=(\$${column}'"\t"')
        fi
    done
    to_print=$(IFS="" ; echo "${to_print[*]}")
    if [[ $file = *.gz ]]; then
	cmd=(gunzip -c $file \| awk "'{print $to_print}'")
    else
        cmd=(awk "'{print $to_print}'" $file)
    fi
    eval "${cmd[@]}" > $output
}


#### CHECK that the base file is of the correct format SNP A1 A2 OR P

for index in ${!base_files[@]};do
    base_file=${base_files[$index]}
    IFS=, read -r -a column <<< "${columns[$index]}"
    pop=${populations[$index]}
    if [[ $base_file = *.gz ]]; then	
	second_row=$(gunzip -c $base_file | head -n 2 | tail -n 1)
    else
	second_row=$(head -n 2 $base_file | tail -n 1)
    fi
    number_of_col=$(echo $second_row| awk '{print NF}')
    [[ ! $number_of_col -eq 5 ]] && [[ -z "${column}" ]] && echo "ERROR: column is for ${base_file} not defined when the base file is not in the correct format or the base file provided do not have 5 columns" && exit 1
    example_id=$(echo $second_row | awk '{print $1}') #make sure that the id is chr:bp

    [[ $number_of_col -eq 5 ]] && [[ ! $example_id = *':'* ]] && echo "The ID of ${base_file} is not in the chr:bp format" && exit 1

    if [[ ! $number_of_col -eq 5 ]] && [[ ! -z ${column} ]]; then
	echo "You want to select 5 columns from the base file ${pop}"
	return_columns $base_file $output_dir/$output_file.$pop.preprocessed_summary_stat ${column[@]}
	[[ ! -f $output_dir/$output_file.$pop.preprocessed_summary_stat ]] && echo $output_dir/$output_file.$pop.preprocessed_summary_stat && exit 1
	base_files[$index]=$output_dir/$output_file.$pop.preprocessed_summary_stat
    fi
done


function join_by { local IFS="$1"; shift; echo "$*"; } #usage joinby [delimiter] [array] => a,b,c

echo "############ Performing PRScsx ##############"


python ~/Desktop/dHCP_genetics/codes/gene_set/PRScsx/PRScsx.py \
        --ref_dir=$ld_folder \
	--bim_prefix=$target_file \
	--sst_file=$(join_by , ${base_files[@]}) \
	--n_gwas=$n_gwas \
	--pop=$(join_by , ${populations[@]}) \
	--phi=1e-2 \
	--meta=TRUE \
	--out_dir=$output_dir \
	--out_name=$output_file

#concatenate the individual files together


for pop in ${populations[@]}; do
    echo concatenate $pop posterior effect sizes
    cat $output_dir/${output_file}_${pop}_pst_eff_a1*.txt > $output_dir/${output_file}_${pop}_all.txt

    awk '{print $1"\t"$1":"$3"\t"$3"\t"$4"\t"$5"\t"$6}' $output_dir/${output_file}_${pop}_all.txt > $output_dir/${output_file}_${pop}_preprocessed.txt

done

#######Calculate the scores######
#####using plink score ####

for pop in ${populations[@]}; do

echo "Calculate PRS score using PLINK for ${pop} population"

plink --bfile $target_file --score $output_dir/${output_file}_${pop}_preprocessed.txt 2 4 6 sum --out $output_dir/${output_file}_${pop}.PRScsx

done
