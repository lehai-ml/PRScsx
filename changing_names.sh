#!/bin/bash

new_output=~/Desktop/dHCP_genetics/dataset/preprocessed_dataset/batch2_HAI/EUR_SAS_EAS/PRS/PScsx_new/
old_output=~/Desktop/dHCP_genetics/dataset/preprocessed_dataset/batch2_HAI/EUR_SAS_EAS/PRS/PScsx/

root=combined_cohort_eur_sas_eas

new_root=batch2_EUR_SAS_EAS_genotyped.50.SCZ

cd $old_output
all_files=$(ls .)

for file in  ${all_files[@]}; do
    file1=$old_output/$file
    file2=$new_output/$new_root${file#$root}

    cmp --silent $file1 $file2 || echo "$file1 and $file2 are different"

done


