awk -F'\t' '$3 == "transcript" { 
    match($0, /gene_id "([^"]+)";/, gene_id); 
    match($0, /transcript_id "([^"]+)";/, transcript_id); 
    print "\"" gene_id[1] "\"; \"" transcript_id[1] "\""; 
}' gtf_file.gtf > extracted_ids.txt
