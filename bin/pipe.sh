declare -a buckets=(100 1000)
declare -a typs=(original normPerLen normPerSeq normPerLenSeq normByLen normBySeq normByLenSeq)

ori=/nfs/users2/cn/sjin/projects/homoplasy/nf_homoplasty
cd $ori

# =============================================
# SCRIPTS


printSh=bin/organize_data/print.sh
raw2deltaPy=bin/organize_data/raw2delta.py
csv2excelPy=bin/organize_data/csv2excel.py
raw2norm2excelcsvPy=bin/organize_data/raw2norm2excelcsv.py

mrdeltaSh=bin/cedric/compute_mrdelta_data.sh
mrdeltaPl=bin/cedric/csv2analysis.pl
mrdeltaPlotGlobalR=bin/cedric/plot_mrdelta_global.R



# =============================================
# RAW & DELTA data: csv & xlsx

# Raw & delta csv
raw_delta_csv () {
echo "# =============================================================="
echo "# Printing RAW & DELTA data"
echo
for bucket in ${buckets[@]}; do
  [[ ! -d data/raw/bucket${bucket} ]] && mkdir -p data/raw/bucket${bucket}
  [[ ! -d data/delta/bucket${bucket} ]] && mkdir -p data/delta/bucket${bucket}
  echo "+print raw bucket${bucket}"
  bash $printSh $bucket
  echo "+print delta bucket${bucket}"
  for i in $(echo "tc sp homo whomo whomo2 ngap ngap2 col len"); do
      file=data/raw/bucket${bucket}/${i}.csv
      filename=${file##*/}
      echo -e "\t+print delta ${filename%.csv} bucket${bucket}"
      python $raw2deltaPy $file data/delta/bucket${bucket}/${filename}
  done
done 
echo 
echo
}

# Raw & delta original & normalized csv & excels
csv2excel_original_norm () {
echo "# =============================================================="
echo "# Manage data (original & normalized) & transform csv to xlsx"
echo 
echo "+norm ${typs[@]}"
for bucket in ${buckets[@]}; do
  [[ ! -d excels/raw/bucket${bucket} ]] && mkdir -p excels/raw/bucket${bucket}
  [[ ! -d excels/delta/bucket${bucket} ]] && mkdir -p excels/delta/bucket${bucket}
  echo "+norm csv2xlsx raw bucket${bucket}"
  python $raw2norm2excelcsvPy $bucket raw
  echo "+norm csv2xlsx delta bucket${bucket}"
  python $raw2norm2excelcsvPy $bucket delta
done
echo
echo
}


# =============================================
# MRDELTA ANALYSIS

mrdelta_analysis () {
# Compute cedric's script and get .tsv files
echo "# =============================================================="
echo "# Compute mrdelta script"
echo
echo "+mrdelta ${typ[@]}"
[[ ! -d data/mrdelta ]] && mkdir data/mrdelta
for bucket in ${buckets[@]}; do
    echo "+mrdelta bucket${bucket}"
    bash $mrdeltaSh $bucket 
done
echo
echo
}

mrdelta_plot () {
# Plot 
# Nratio or deltaacc vs mrdelta
echo "# =============================================================="
echo "# Plot mrdelta"
echo
echo "+plot mrdelta global ${typ[@]}"
[[ ! -d plots/mrdelta/global/bucket100_1000 ]] && mkdir -p plots/mrdelta/global/bucket100_1000
Rscript $mrdeltaPlotGlobalR
}


#raw_delta_csv
#csv2excel_original_norm
#mrdelta_analysis
mrdelta_plot
