table_input=$1
blastn_database=$2

test -e ${blastn_database} && rm ${blastn_database}
cat ${table_input} |while read database lable
do
  echo $database $lable
  sed "s/>/>$lable:/" $database >>${blastn_database}
done

makeblastdb  -in ${blastn_database} -dbtype nucl 
