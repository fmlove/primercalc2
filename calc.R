# Sanger seq primer calculator

target_length = 20
length_var = 2
GC_content_goal = 50
threeprime_GC_goal = TRUE
melting_temp_goal = 55

#TODO - handle split without separate sequence
complement <- function(sequence = NULL, split = NULL){
  if(missing(sequence) & missing(split)){stop("A DNA sequence must be provided.")}
  if(missing(split)){split = strsplit(sequence, '')}
  comp = sapply(split, function(b){
    if(b == 'A'){'T'}
    else if(b == 'T'){'A'}
    else if(b == 'C'){'G'}
    else if(b == 'G'){'C'}
  })
  return(paste0(comp))
}

antiparallel <- function(sequence, split = NULL){
  if(missing(sequence) & missing(split)){stop("A DNA sequence must be provided.")}
  if(missing(split)){split = strsplit(sequence, '')}
  comp = complement(split = split)
  ap = rev(comp)
  return(paste0(ap))
}



check_DNA <- function(object){#TODO - allow ambiguous bases?
  
}

setClass('DNA', 
         representation(sequence = 'character', 
                        complement = 'character', 
                        antiparallel = 'character', 
                        length = 'numeric', 
                        GC_percent = 'numeric', 
                        melting_temp = 'numeric',
                        composition = 'table'),
         validity = check_DNA
         )

DNA <- function(sequence){
  sequence = toupper(sequence)
  split = unlist(strsplit(sequence, ''))
  complement = complement(split = split)
  antiparallel = antiparallel(split = split)
  length = nchar(sequence)
  seq.table = table(split)
  GC_percent = ((seq.table['C'] + seq.table['G'])/length)*100
  melting_temp = seq.table['A']*2 + seq.table['T']*2 + seq.table['G']*4 + seq.table['C']*4
  
  new('DNA', sequence = sequence, 
      complement = complement,
      antiparallel = antiparallel,
      length = length,
      GC_percent = GC_percent,
      melting_temp = melting_temp,
      composition = seq.table
      )
}

calculate_seq_primer <- function(sequence, rev = FALSE){
  seq = DNA(sequence)@sequence
  if(rev){
    seq@antiparallel
  }
  candidates = list()
  start = 1
  len = target_length - length_var
  while(len <= target_length + length_var){
    while(start + len <= nchar(seq)){
      s = substr(seq, start, start + (len - 1))
      candidates = c(candidates, DNA(s))
      start = start + 1
    }
    start = 1
    len = len + 1
  }
  data.frame(id = 1:length(candidates),
             s = sapply(candidates, function(c){c@sequence}),
             GC_percent = sapply(candidates, function(c){c@GC_percent}),
             mt = sapply(candidates, function(c){c@melting_temp}),
             GC_end = grepl("(G|C)$", sapply(candidates, function(c){c@sequence}))
             )
}

