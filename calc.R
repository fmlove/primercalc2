# Sanger seq primer calculator

target_length = 20
length_var = 2
GC_content_goal = 52
GC_content_var = 2
threeprime_GC_goal = TRUE
melting_temp_goal = 52
melting_temp_var = 2
hairpin_min = 4
  
#TODO - handle split without separate sequence
complement <- function(sequence = NULL, split = NULL){
  if(missing(sequence) & missing(split)){stop("A DNA sequence must be provided.")}
  if(missing(split)){split = unlist(strsplit(sequence, ''))}
  comp = sapply(split, function(b){
    if(b == 'A'){'T'}
    else if(b == 'T'){'A'}
    else if(b == 'C'){'G'}
    else if(b == 'G'){'C'}
  })
  return(paste0(comp, collapse = ''))
}

antiparallel <- function(sequence, split = NULL){
  if(missing(sequence) & missing(split)){stop("A DNA sequence must be provided.")}
  if(missing(split)){split = unlist(strsplit(sequence, ''))}
  comp = complement(split = split)
  ap = rev(unlist(strsplit(comp, '')))
  return(paste0(ap, collapse = ''))
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
  GC_percent = (sum(seq.table['C'], seq.table['G'], na.rm = TRUE)/length)*100
  melting_temp = sum(seq.table['A']*2, seq.table['T']*2, seq.table['G']*4, seq.table['C']*4, na.rm = TRUE)
  
  new('DNA', sequence = sequence, 
      complement = complement,
      antiparallel = antiparallel,
      length = length,
      GC_percent = GC_percent,
      melting_temp = melting_temp,
      composition = seq.table
      )
}


check_self_annealing <- function(sequence = NULL, split = NULL){
  if(missing(sequence) & missing(split)){stop("A DNA sequence must be provided.")}
  if(missing(split)){split = unlist(strsplit(sequence, ''))}
  split = toupper(split)
  seq = paste0(split, collapse = '')
  start = 1
  hairpin = FALSE
  starts = list()
  while(start + (hairpin_min - 1) + hairpin_min <= length(split)){
    scan = antiparallel(split = split[start:(start + (hairpin_min - 1))])
    s = substr(seq, start + hairpin_min, nchar(seq))
    if(grepl(scan, s)){hairpin = TRUE; starts = c(starts, start)}
    start = start + 1
  }
  return(hairpin)
}


calculate_seq_primer <- function(sequence, rev = FALSE){
  seq.DNA = DNA(sequence)
  seq = seq.DNA@sequence
  if(rev){
    seq = seq.DNA@antiparallel
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
  
  cdf = data.frame(id = 1:length(candidates),
             s = sapply(candidates, function(c){c@sequence}),
             length = sapply(candidates, function(c){nchar(c@sequence)}),
             GC_percent = sapply(candidates, function(c){c@GC_percent}),
             mt = sapply(candidates, function(c){c@melting_temp}),
             GC_end = grepl("(G|C)$", sapply(candidates, function(c){c@sequence})),
             stringsAsFactors = FALSE
  )
  
  mt_in_range = sapply(candidates, function(c){c@melting_temp >= melting_temp_goal - melting_temp_var & c@melting_temp <= melting_temp_goal + melting_temp_var})
  gc_in_range = sapply(candidates, function(c){c@GC_percent >= GC_content_goal - GC_content_var & c@GC_percent <= GC_content_goal + GC_content_var})
  gc_end = grepl("(G|C)$", sapply(candidates, function(c){c@sequence}))
  cdf$no_dimers = sapply(candidates, function(c){!check_self_annealing(c@sequence)})
  
  cdf$good_candidate = mt_in_range & gc_in_range & gc_end & cdf$no_dimers
  
  return(cdf)
}

