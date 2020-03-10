

declare -a aligner=(CLUSTALO MAFFT-FFTNS1 MAFFT-GINSI MAFFT-SPARSECORE)

declare -f folders=(results_fullTree_CO_bucket100 results_fullTree_FFTNS1_bucket100 results_fullTree_GINSI_bucket100 results_fullTree_SPARSE_bucket100)

declare -a familyName=(ltn     il8     az      kringle cryst   DEATH   cah     mmp     rub     ghf10   tgfb    sodcu   KAS     DMRL_synthase   tms     GEL     kunitz  Sulfotransfer   mofe    Ald_Xan_dh_2    ghf5    phc     aadh    annexin serpin  cytb    asp     oxidored_q6     hpr     hormone_rec     hr      tim     glob    ace     cys     ghf1    sodfe   peroxidase      uce     flav    HMG_box OTCace  msb     icd     proteasome      cyclo   LIM     HLH     ldh     subt    int     lyase_1 gpdh    egf     blm     gluts   myb_DNA-binding tRNA-synt_2b    biotin_lipoyl   hom     ghf13   aldosered       hla     Rhodanese       PDZ     blmb    rhv     p450    adh     aat     rrm     Acetyltransf    sdr     zf-CCHH rvp)
##declare -a familyName=(seatoxin)

declare -a tree=(codnd dpparttreednd1 dpparttreednd2 dpparttreednd2size fastaparttreednd fftns1dnd fftns1dndmem fftns2dnd fftns2dndmem mafftdnd parttreednd0 parttreednd1 parttreednd2 parttreednd2size FAMSA CLUSTALO-RANDOM)
##declare -a tree=(codnd dpparttreednd1 dpparttreednd2)

declare -a tcs=(tc sp col)
declare -a homos=(homo w_homo w_homo2 len ngap ngap2)


#familylist=/users/cn/sjin/projects/homoplasy/homfam_info/list_of_homfam_families

printf "family aligner tree1 tree2 tc sp col homo w_homo w_homo2 len ngap ngap2\n"

#for fam in $(cat $familylist); do
for fam in ${familyName[@]}; do
    i=0
    for aln in ${aligner[@]}; do
        fold=${folders[$i]}
        ntr=0
        for tr1 in ${tree[@]}; do
        ntr=$((ntr+1))
        for tr2 in ${tree[@]:$ntr}; do
            printf "$fam $aln $tr1 $tr2 "
            for tc in ${tcs[@]}; do
                fil1=${fold}/individual_scores/${fam}.dpa_align.100.${aln}.${tr1}.${tc}
                fil2=${fold}/individual_scores/${fam}.dpa_align.100.${aln}.${tr2}.${tc}
                if [[ -s $fil1 && -s $fil2 ]]; then
                    val1=$(awk 'NR==1{printf $1" "}' ${fil1})
                    val2=$(awk 'NR==1{printf $1" "}' ${fil2})
                    python -c "print $val1 - $val2" | tr "\n" " "
                else
                    printf "NA "
                fi
            done
            for hom in ${homos[@]}; do
                fil1=${fold}/alignments/${fam}.dpa_100.${aln}.with.${tr1}.tree.${hom}
                fil2=${fold}/alignments/${fam}.dpa_100.${aln}.with.${tr2}.tree.${hom}
                if [[ -s $fil1 && -s $fil2 ]]; then
                    val1=$(awk 'NR==1{printf $1" "}' ${fil1} | tr -d "l")
                    val2=$(awk 'NR==1{printf $1" "}' ${fil2} | tr -d "l")
                    if [[ $hom == "ngap2" ]]; then 
                        python -c "print $val1 - $val2" 
                    else
                        python -c "print $val1 - $val2" | tr "\n" " "
                    fi
                else
                    if [[ $hom == "ngap2" ]]; then 
                        printf "NA\n"
                    else
                        printf "NA "
                    fi
                fi
            done
        done
        done
        i=$((i+1))
    done
done