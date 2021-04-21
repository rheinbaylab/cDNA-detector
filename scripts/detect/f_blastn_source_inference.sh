input_file=$1
db=$2
blastn -query ${input_file} -db ${db} -task "blastn-short" -word_size 6 -max_target_seqs 1 -outfmt 7 -ungapped  > ${input_file}.blastn
grep "^# Fields" ${input_file}.blastn  |sed 's/# Fields: //;s/, /\t/g;s/%/pct/' |sed 's/ /_/g' |head -n 1 >${input_file}.blastn.tsv
grep -v "^#" ${input_file}.blastn >>${input_file}.blastn.tsv
rm ${input_file}.blastn


