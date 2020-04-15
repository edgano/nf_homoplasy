t_coffee -reg -reg_method dynamic_msa \
         -seq ${seqs} \
         -reg_nseq ${bucket_size} \
         -reg_dynamic 10000 \
         -reg_homoplasy \
         -blast_server LOCAL \
         -protein_db=/Users/edgargarriga/CBCRG/bin/uniref50.fasta \
         -outfile ${id}.reg_align.${bucket_size}.${align_method}.with.${tree_method}.tree.aln
