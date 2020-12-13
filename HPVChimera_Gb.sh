blast(){
echo "----------blast--------"
filename=$2
folder=$1
hpv_database=$3

blastn -db $hpv_database -query "$workspace"/"$folder"/"$filename"  -task blastn -evalue 0.01 -max_target_seqs 10 -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore" -out "$workspace"/"$folder"/"$filename".out
cat "$workspace"/"$folder"/"$filename".out | sort -nrk 13 >"$workspace"/"$folder"/"$filename".ordered
}


compare_similarity_hpv(){

folder=$1
file=$2

if [ `ls -l "$workspace"/"$folder"/"$file" | awk '{print $5}'` -eq 0 ]
then
   filename=$(echo $2| awk -F ".ordered" '{print $1}')
   echo "File: $file -> Chimera: Contig sequence does not show any similarity to any HPV type given: Blaclisted"
   echo "$file -> Chimera: Contig sequence does not show any similarity to any HPV type given">>$result_file
   echo "$file" "Step1" >>$blacklisted
else
  bit_score=`cat "$workspace"/"$folder"/"$file" | head -1 | awk '{print $3}'`
  if (( $(echo "$bit_score < 85.0" |bc -l) ))
  then
    filename=$(echo $2| awk -F ".ordered" '{print $1}')
    echo "File: $file -> Chimera: Contig sequence does not show >85% sequence identity to any HPV type given.: Blacklisted"
    echo "$file -> Chimera: Contig sequence does not show >85% sequence identity to any HPV type given.">>$result_file
    echo "$file" "Step2">>$blacklisted
  fi
fi


}


compare_similarity_hpv_div(){

folder=$1
file=$2

if [ `ls -l "$workspace"/"$folder"/"$file" | awk '{print $5}'` -eq 0 ]
then
   filename=$(echo $2| awk -F ".ordered" '{print $1}')
   echo "File: $file -> Chimera: At least one contig segment sequence does not show any similarity to any HPV type given: Blaclisted"
   echo "$file -> Chimera: At least one contig segment does not show any similarity to any HPV type given">>$result_file
   echo "$file" "Step 8" >>$blacklisted
else
  bit_score=`cat "$workspace"/"$folder"/"$file" | head -1 | awk '{print $3}'`
  if (( $(echo "$bit_score < 85.0" |bc -l) ))
  then
    filename=$(echo $2| awk -F ".ordered" '{print $1}')
    echo "File: $file -> Chimera: At least one contig segment does not show >85% sequence identity to any HPV type given.: Blacklisted"
    echo "$file -> Chimera: At least one contig segment does not show >85% sequence identity to any HPV type given.">>$result_file
    echo "$file" "Step 9">>$blacklisted
  fi
fi


}





check_coverage()
{
  filename="$(basename $folder)"
  isblacklisted=`cat "$blacklisted" | grep "$filename"| wc -l`
  echo "----check coverage-----"
  filename=$2
  folder=$1
  filename=$(echo $2| awk -F ".ordered" '{print $1}')
  filename="$(basename $folder)"
  isblacklisted=`cat "$blacklisted" | grep "$filename"| wc -l`
  if [ $isblacklisted -eq "1" ]; then
    echo "$filename is blacklisted, this file will not continue"
  else
    filename_unique=$2"_unique.out"
    if [ -s "$workspace"/"$folder"/"$filename_unique" ]
    then
      value=$(cat "$workspace"/"$folder"/"$filename_unique" | awk -F " " {'print $4'})
      echo "----------- $extension has some data. Value: $value"
      echo "$value" "-lt" "$value_to_discriminate"  $filename
      if [ "$value" -lt "$value_to_discriminate" ]; then
        echo "$filename -> Chimera: $value: Blacklisted"
        echo "$filename -> Chimera: The top hit alignment does not cover >60% of contigs sequence ">>$result_file
        echo "$filename" "Step 3">>$blacklisted
      else
        echo "$filename_unique -> NO Chimera: $value"
      fi
   else
     echo "$filename_unique -> Chimera: The top hit alignment does not cover >60% of contigs sequence"
     echo "$filename_unique -> Chimera: The top hit alignment does not cover >60% of contigs sequence">>$result_file
     echo "$filename" "Step 4">>$blacklisted
   fi
 fi
}

check_coverage_div()
{
  filename="$(basename $folder)"
  isblacklisted=`cat "$blacklisted" | grep "$filename"| wc -l`
  echo "----check coverage-----"
  filename=$2
  folder=$1
  filename=$(echo $2| awk -F ".ordered" '{print $1}')
  filename="$(basename $folder)"
  isblacklisted=`cat "$blacklisted" | grep "$filename"| wc -l`
  if [ $isblacklisted -eq "1" ]; then
    echo "$filename is blacklisted, this file will not continue"
  else
    filename_unique=$2"_unique.out"
    if [ -s "$workspace"/"$folder"/"$filename_unique" ]
    then
      value=$(cat "$workspace"/"$folder"/"$filename_unique" | awk -F " " {'print $4'})
      echo "----------- $extension has some data. Value: $value"
      echo "$value" "-lt" "$value_to_discriminate_div"  $filename
      if [ "$value" -lt "$value_to_discriminate_div" ]; then
        echo "$filename -> Chimera: $value: Blacklisted"
        echo "$filename -> Chimera: The top hit alignment does not cover >70% of contigs segment sequence ">>$result_file
        echo "$filename" "Step 10">>$blacklisted
      else
        echo "$filename_unique -> NO Chimera: $value"
      fi
   else
     echo "$filename_unique -> Chimera: The top hit alignment does not cover >70% of contigs segment sequence"
     echo "$filename_unique -> Chimera: The top hit alignment does not cover >70% of contigs segment sequence">>$result_file
     echo "$filename" "Step 11">>$blacklisted
   fi
 fi
}


get_only_first_line_unique()
{
echo "--get_only_first_line_unique-----------"
filename=$2
folder=$1
filename_ordered=$2".ordered"
filename="$(basename $folder)"
isblacklisted=`cat "$blacklisted" | grep "$filename"| wc -l`
if [ $isblacklisted -eq "1" ]; then
  echo "$filename is blacklisted, this file will not continue"
else
  echo "$filename is not blacklisted, we can continue"
   echo  ""$workspace"/"$folder"/"$filename".fa.ordered"
   bit_score_first=$(cat "$workspace"/"$folder"/"$filename".fa.ordered | head -1 | awk -F " " '{print $13}')
   HPV_type_first=$(cat "$workspace"/"$folder"/"$filename".fa.ordered | head -1 | awk -F " " '{print $2}')
   count=1
   echo $count>$workspace/count_txt.log
   echo "--------------------------"
   echo "Getting first value in order to compare BitScore"
   cat "$workspace"/"$folder"/"$filename.fa.ordered"|tail -n +2 | while read line;do
     bit_score=$(echo $line | awk -F " " '{print $13}')
     HPV_type=$(echo $line | awk -F " " '{print $2}')
     echo $bit_score  $bit_score_first
     if [ "$bit_score" == "$bit_score_first" ]; then
       if [ "$HPV_type" == "$HPV_type_first" ]; then
         echo "HPV type scored with same bit"
         count=1
       else
         count=$((count+1))
         echo $count>$workspace/count_txt.log
       fi
     else
       count=1
       break
     fi
     echo $count>$workspace/count_txt.log
   done
   count=$(cat "$workspace"/count_txt.log)
   rm -rf $workspace/count_txt.log
   if [ $count -eq "1" ]; then
     echo -e "Number of Bitscore found" $count
     cat ""$workspace"/"$folder"/"$filename.fa.ordered"" | head -$count >"$workspace"/"$folder"/"$filename"_unique.out
     echo -e "Unique head for $filename"
   else
     echo "$filename -> Chimera: The sequence aligns to 2 different HPV types with same bitscore">>$result_file
     echo "$filename -> Chimera: The sequence aligns to 2 different HPV types with same bitscore: Blacklisted"
     echo "$filename" "Step 5">>$blacklisted
   fi
fi


}

get_only_first_line_unique_div()
{
echo "--get_only_first_line_unique_div-----------"
echo $file
file=$2
folder=$1
filename_ordered=$2".ordered"
filename="$(basename $folder)"
isblacklisted=`cat "$blacklisted" | grep "$filename"| wc -l`
if [ $isblacklisted -eq "1" ]; then
  echo "$filename is blacklisted, this file will not continue"
else
   echo "$filename is not blacklisted, we can continue"
   bit_score_first=$(cat "$workspace"/"$folder"/"$file" | head -1 | awk -F " " '{print $13}')
   HPV_type_first=$(cat "$workspace"/"$folder"/"$file" | head -1 | awk -F " " '{print $2}')
   echo $bit_score_first
   count=1
   echo $count>$workspace/count_txt.log
   echo "--------------------------"
   echo "Getting first value in order to compare BitScore"
   echo "$workspace"/"$folder"/"$file"
   cat "$workspace"/"$folder"/"$file"|tail -n +2 | while read line;do
     bit_score=$(echo $line | awk -F " " '{print $13}')
     HPV_type=$(echo $line | awk -F " " '{print $2}')
     echo $bit_score  $bit_score_first
     if [ "$bit_score" == "$bit_score_first" ]; then
       if [ "$HPV_type" == "$HPV_type_first" ]; then
         count=1
         echo "HPV type scored with same bit"
       else
         count=$((count+1))
         echo $count>$workspace/count_txt.log
       fi
     else
       count=1
       break
     fi
     echo $count>$workspace/count_txt.log
   done
   count=$(cat "$workspace"/count_txt.log)
   rm -rf $worspace/count_txt.log
   if [ $count -eq "1" ]; then
     echo -e "Number of Bitscore found" $count
     cat "$workspace"/"$folder"/"$file" | head -$count > "$workspace"/"$folder"/"$file"_unique.out
     echo -e "Unique head for $filename"
   else
     echo "$filename -> Chimera: At least one sequence segment aligns to 2 different HPV types with same bitscore">>$result_file
     echo "$filename -> Chimera: At least one sequence segment aligns to 2 different HPV types with same bitscore: Blacklisted"
     echo "$filename" "Step 12">>$blacklisted
   fi
fi
}




divide_in_three(){
echo "----------divide_in_three------------"

file=$2
folder=$1
isblacklisted=`cat "$blacklisted" | grep "$filename"| wc -l`
if [ $isblacklisted -eq "1" ]; then
  echo "$filename is blacklisted, this file will not continue"
else
  echo "$filename is not blacklisted, we can continue"
  extension="${filename%.*}"
  echo -e "$filename is a candidate"
  long=$(cat "$workspace"/"$folder"/"$filename".fa | tail -n +2|sed "/\w+/g" |wc -m)
  echo $long
  sequence=$(cat "$workspace"/"$folder"/"$filename".fa | tail -n +2)
  L1=$(expr $long / 3)
  L11=$(expr $L1 - 1)
  L2=$(expr $long \* 2)
  L2=$(expr $L2 / 3)
  echo $filename --- $long "----" $L1 "----" $L2
  var1="${sequence:0:$L1}"
  var2="${sequence:$L1:$L2-$L1}"
  var3="${sequence:$L2}"
  echo -e ">"$extension"_1">"$workspace"/"$folder"/"$file"_1.fasta
  echo -e $var1>>"$workspace"/"$folder"/"$file"_1.fasta
  echo -e ">"$extension"_2">"$workspace"/"$folder"/"$file"_2.fasta
  echo -e $var2>>"$workspace"/"$folder"/"$file"_2.fasta
  echo -e ">"$extension"_3">"$workspace"/"$folder"/"$file"_3.fasta
  echo -e $var3>>"$workspace"/"$folder"/"$file"_3.fasta
fi
}



blast_85_div()
{

echo "-----------blast_85_div------"
file=$2
folder=$1

isblacklisted=`cat "$blacklisted" | grep "$filename"| wc -l`

if [ $isblacklisted -eq "1" ]; then
  echo "$filename is blacklisted, this file will not continue"
else
  echo -e "blastn for $filename"
  #blastn -db $hpv_database -query "$workspace"/"$folder"/"$file"  -task blastn -evalue 0.01 -max_target_seqs 10 -word_size 8 -dust no -soft_masking false -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore" -out "$workspace"/"$folder"/"$file""_concoverageb.out"

  HPV_type_top=$(cat "$workspace"/"$folder"/"$folder"_unique.out | awk -F " " '{print $2}')
  HPV_type_top_file_path="$folder_database_to_single_fasta"/"$HPV_type_top".fa
  echo $HPV_type_top_file_path
  blastn  -query "$workspace"/"$folder"/"$file"  -subject "$HPV_type_top_file_path" -task blastn  -max_target_seqs 10 -word_size 8 -evalue 0.01  -dust no -soft_masking false -outfmt "6 qseqid sseqid pident qcovs length mismatch gapopen qstart qend sstart send evalue bitscore" -out "$workspace"/"$folder"/"$file""_concoverageb.out"
  cat "$workspace"/"$folder"/"$file""_concoverageb.out" | sort -nrk 13 >"$workspace"/"$folder"/"$file""_concoverageb.out".ordered
fi

}


last_check()
{
echo "----------last check---------"
isblacklisted=`cat "$blacklisted" | grep "$filename"| wc -l`
folder=$1

if [ $isblacklisted -eq "1" ]; then
  echo "$filename is blacklisted, this file will not continue"
else

 echo "$folder -> $HPV_type_original, no chimera detected">>$result_file
 echo "$folder -> $HPV_type_original, no chimera detected"
fi
}




extract_id_indivial_files() {
cat $fasta_file | sed 's/ //g' |sed 's/;//g' |sed 's/{//g'|sed 's/}//g'|sed 's/\]//g'|sed 's/_//g'|sed 's/\[//g'|sed 's/://g'|sed 's/\,//g' |sed 's/(//g'| sed 's/)//g'|sed 's/\.//g' | sed 's/\-//g' | sed 's/\|//g' |   sed 's/\///g' | sed 's/|//g'|cut -c -100 |sed $'s/[^[:print:]\t]//g' |awk '{if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}print $0 > filename}'

cd  $folder_database_to_single_fasta
cat $folder_database_to_single_fasta/$hpv_database_file_name| sed 's/ //g' |sed 's/;//g' |sed 's/{//g'|sed 's/}//g'|sed 's/\]//g'|sed 's/_//g'|sed 's/\[//g'|sed 's/://g'|sed 's/\,//g' |sed 's/(//g'| sed 's/)//g'|sed 's/\.//g' | sed 's/\-//g' | sed 's/\|//g' |   sed 's/\///g' | sed 's/|//g'|cut -c -100 |sed $'s/[^[:print:]\t]//g' |awk '{if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")}print $0 > filename}'

cd $folder_to_start/chimera/workspace
rm -rf $folder_database_to_single_fasta/$hpv_database_file_name
echo -e "Fasta file well divided"
}

clean_gargabe()
{
echo "Clean Gargabe"
find  "$workspace"/"$folder"/ | grep -v unique | xargs rm -rf
}

if_fasta_file() {
  file=$1
  extension="${file##*.}"
  if [ ! -f "$file" ] || [ $extension != "fasta" ]; then
    echo "HPV Chimera check failed."
    echo "This is not a .fasta file or the $file does exist."
    exit $?
  fi
}

checkExit()
{
 if [ "$?" -ne 0 ]; then
    set -e
    echo Command $1 exited with abnormal status
    echo "HPV Chimera check failed."
    exit $?;
  fi
}


help()
{
  echo "***********************************"
  echo -e "Chimera Check"
  echo -e "Three arguemtns: "
  echo -e "1. Assembled contigs to analyze in fasta file. Full Path."
  echo -e "2. Folder to work. Full Path."
  echo -e "3. HPV database fasta file. Full Path."
  echo -e "''''''''''''''''''"
}
header()
{
  echo "***********************************"
  echo -e "HPV Chimera Check"
  echo "************************************"

}
final()
{
  echo "HPV Chimera Check Finished"
  echo "-------------------------"
  echo "Check result in: "$folder_to_start"chimera/HPVchimera_results.txt"
  echo "Check files blacklisted in "$blacklisted""
}

if [ "$#" -ne 3 ]; then
  echo "Illegal number of parameters"
  echo "HPV Chimera check failed."
  help
  exit 0
fi



if [ "$1" == "-h" ] ; then
  echo "Usage: `basename $0` [-h]"
  help
  exit 0
fi

fasta_file=$1
fasta_file_name="$(basename $fasta_file)"
folder_to_start=$2
hpv_database=$3

set -e
mkdir -p $folder_to_start/chimera
val_identity_ninety=90.000
val_identity_five=5.000
root_folder=$folder_to_start/chimera
workspace=$root_folder/workspace
value_to_discriminate=60
value_to_discriminate_div=70
result_file=$folder_to_start/chimera/HPVchimera_results.txt
rm -rf $result_file
touch $result_file
blacklisted=$folder_to_start/chimera/HPVchimera_blacklisted.txt
rm -rf $blacklisted
touch $blacklisted
if_fasta_file $fasta_file
if_fasta_file $hpv_database
cp $fasta_file $root_folder
cd $folder_to_start/chimera/

rm -rf workspace
mkdir -p workspace
mkdir -p database_to_single_fasta
cp $hpv_database $root_folder/database_to_single_fasta/
hpv_database_file_name="$(basename -- $hpv_database)"
folder_database_to_single_fasta=$root_folder/database_to_single_fasta/
cp $fasta_file workspace
cd workspace
header
extract_id_indivial_files
rm -rf $folder_to_start/chimera/workspace/$fasta_file_name


#cat $workspace/fasta_file_loop.fa | grep gi | sed 's/^.//' | while read filename;do
ls | grep -v fasta| while read filename;do
  folder_to_create="${filename%.*}"
  mkdir -p $folder_to_create
  mv $filename "$folder_to_create"/"$filename"
  blast "$folder_to_create" "$filename" "$hpv_database"
  compare_similarity_hpv "$folder_to_create" "$filename".ordered
  get_only_first_line_unique "$folder_to_create" "$filename"
  check_coverage "$folder_to_create" "$filename"
  divide_in_three "$folder_to_create" "$filename"
  blast_85_div "$folder_to_create" "$filename"_1.fasta
  blast_85_div "$folder_to_create" "$filename"_2.fasta
  blast_85_div "$folder_to_create" "$filename"_3.fasta
  compare_similarity_hpv_div "$folder_to_create" "$filename"_1.fasta_concoverageb.out.ordered
  compare_similarity_hpv_div "$folder_to_create" "$filename"_2.fasta_concoverageb.out.ordered
  compare_similarity_hpv_div "$folder_to_create" "$filename"_3.fasta_concoverageb.out.ordered
  get_only_first_line_unique_div "$folder_to_create" "$filename"_1.fasta_concoverageb.out.ordered
  get_only_first_line_unique_div "$folder_to_create" "$filename"_2.fasta_concoverageb.out.ordered
  get_only_first_line_unique_div "$folder_to_create" "$filename"_3.fasta_concoverageb.out.ordered
  check_coverage_div "$folder_to_create" "$filename"_1.fasta_concoverageb.out.ordered
  check_coverage_div "$folder_to_create" "$filename"_2.fasta_concoverageb.out.ordered
  check_coverage_div "$folder_to_create" "$filename"_3.fasta_concoverageb.out.ordered
  last_check "$folder_to_create" "$filename"
done
final
