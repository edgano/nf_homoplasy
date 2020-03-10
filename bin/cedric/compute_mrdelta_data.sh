bucket=bucket"$1"
declare -a mrdeltas=(0 0.05 0.1 0.15 0.2 0.25 0.3 0.35 0.4 0.45 0.5 0.6 0.7 0.8 0.9 1)
declare -a scores=(tc sp)
declare -a metrics=(homo whomo whomo2 ngap ngap2 homoPerLen whomoPerLen whomo2PerLen ngapPerLen ngap2PerLen homoByLen whomoByLen whomo2ByLen ngapByLen ngap2ByLen homoPerSeq whomoPerSeq whomo2PerSeq ngapPerSeq ngap2PerSeq homoBySeq whomoBySeq whomo2BySeq ngapBySeq ngap2BySeq homoPerLenSeq whomoPerLenSeq whomo2PerLenSeq ngapPerLenSeq ngap2PerLenSeq homoByLenSeq whomoByLenSeq whomo2ByLenSeq ngapByLenSeq ngap2ByLenSeq)
#declare -a metrics=(whomo whomo2 whomoPerLen whomo2PerLen whomoByLen whomo2ByLen whomoPerSeq whomo2PerSeq whomoBySeq whomo2BySeq whomoPerLenSeq whomo2PerLenSeq whomoByLenSeq whomo2ByLenSeq)
declare -a aligners=(CLUSTALO MAFFT-FFTNS1 MAFFT-GINSI MAFFT-SPARSECORE)
#declare -a familyName=(ltn     il8     az      kringle cryst   DEATH   cah     mmp     rub     ghf10   tgfb    sodcu   KAS     DMRL_synthase   tms     GEL     kunitz  Sulfotransfer   mofe    Ald_Xan_dh_2    ghf5    phc     aadh    annexin serpin  cytb    asp     oxidored_q6     hpr     hormone_rec     hr      tim     glob    ace     cys     ghf1    sodfe   peroxidase      uce     flav    HMG_box OTCace  msb     icd     proteasome      cyclo   LIM     HLH     ldh     subt    int     lyase_1 gpdh    egf     blm     gluts   myb_DNA-binding tRNA-synt_2b    biotin_lipoyl   hom     ghf13   aldosered       hla     Rhodanese       PDZ     blmb    rhv     p450    adh     aat     rrm     Acetyltransf    sdr     zf-CCHH rvp)
familyList=$(cat /users/cn/sjin/projects/homoplasy/homfam_info/combinedSeqs)

csv2analyse=bin/cedric/csv2analyse.pl
cedrictxt2colcsv=bin/cedric/cedrictxt2colcsv.py


# Create files
compute_homoplasy () {
    ori=data/raw/${bucket}  # Folder where the raw data files are stored
    outori=data/mrdelta
    for mrdelta in ${mrdeltas[@]}; do
        folder=${outori}/${bucket}/mrdelta${mrdelta}
        [[ ! -d $folder ]] && mkdir -p $folder
        for score in ${scores[@]}; do
            for metric in ${metrics[@]}; do
                echo "perl $csv2analyse -mrdelta $mrdelta -score $score -metrics $metric -dir $ori > ${folder}/${score}_${metric}.txt"
                perl $csv2analyse -mrdelta $mrdelta -score $score -metrics $metric -dir $ori > ${folder}/${score}_${metric}.txt
                python $cedrictxt2colcsv ${folder}/${score}_${metric}.txt
            done
        done
    done
}

# Merge data 
# This will be used to plot
merge_data () {
    echo "+mrdelta merge $bucket"
    ori=data/mrdelta/${bucket}
    outori=data/mrdelta  # Where mrdelta merged output tables are stored
    outfam=${outori}/out.${bucket}.fam.tsv
    outaln=${outori}/out.${bucket}.aln.tsv
    outglobal=${outori}/out.${bucket}.global.tsv
    printf "mrdelta\tscore\tmetric\tfamily\taligner\tNratio\tNnumber\tusedNumber\tunusedRatio\tunusedNumber\ttotal\n" > $outfam
    printf "mrdelta\tscore\tmetric\taligner\tNratio\tunusedRatio\tnfam\tminacc\tmaxacc\tdeltaacc\n" > $outaln
    printf "mrdelta\tscore\tmetric\taligner\tNratio\tunusedRatio\tminacc\tmaxacc\tdeltaacc\n" > $outglobal
    echo -e "\t+mrdelta merge global $bucket"
    for mrdelta in ${mrdeltas[@]}; do
        for score in ${scores[@]}; do
            for metric in ${metrics[@]}; do
                awk -v mrdelta=$mrdelta -v score=$score -v metric=$metric 'NR!=1{print mrdelta"\t"score"\t"metric"\tglobal\t"$0}' ${ori}/mrdelta${mrdelta}/${score}_${metric}_global.tsv >> $outglobal 
            done
        done
    done
    echo -e "\t+mrdelta merge aln $bucket"
    for mrdelta in ${mrdeltas[@]}; do
        for score in ${scores[@]}; do
            for aligner in ${aligners[@]}; do
                for metric in ${metrics[@]}; do
                awk -v mrdelta=$mrdelta -v score=$score -v metric=$metric -v aligner=$aligner 'NR!=1&&$1==aligner{print mrdelta"\t"score"\t"metric"\t"$0}' ${ori}/mrdelta${mrdelta}/${score}_${metric}_aln.tsv >> $outaln
                done
            done
        done
    done
    echo -e "\t+mrdelta merge fam $bucket"
    for mrdelta in ${mrdeltas[@]}; do
        for score in ${scores[@]}; do
            for metric in ${metrics[@]}; do
                awk -v mrdelta=$mrdelta -v score=$score -v metric=$metric 'NR!=1{print mrdelta"\t"score"\t"metric"\t"$0}' ${ori}/mrdelta${mrdelta}/${score}_${metric}_fam.tsv >> $outfam
            done
        done
    done
    # for mrdelta in ${mrdeltas[@]}; do
    #     for score in ${scores[@]}; do
    #         #for fam in ${familyName[@]}; do
    #         for fam in $familyList; do
    #         for aligner in ${aligners[@]}; do
    #             for metric in ${metrics[@]}; do
    #                 awk -v mrdelta=$mrdelta -v score=$score -v metric=$metric -v fam=$fam -v aligner=$aligner 'NR!=1&&$1==fam&&$2==aligner{print mrdelta"\t"score"\t"metric"\t"$0}' ${ori}/mrdelta${mrdelta}/${score}_${metric}_fam.tsv >> $outfam
    #             done
    #         done
    #         done
    #     done
    # done
}


compute_homoplasy
merge_data

