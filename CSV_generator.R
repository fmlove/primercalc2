#oligo CSV generation
#manual - import CSV with names and sequences (as 'seq')

#properties
generate_properties <- function(df, fileout){
  seqs = sapply(df$seq, DNA)
  df$Tm = sapply(seqs, function(s){s@melting_temp})
  df$GC = sapply(seqs, function(s){s@GC_percent})
  write.csv(df, fileout)
}


