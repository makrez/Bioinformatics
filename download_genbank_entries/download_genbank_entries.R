extract_sequences_genbank <- function(organism,
                                      sequence_length = NULL,
                                      other_search_words = NULL,
                                      gene = NULL,
                                      filename = NULL,
                                      max_entries = 40,
                                      datatype = "fasta",
                                      database = "nucleotide"
                                      ){
  # script to download sequences form Genbank using the package rentrez.
  # Args:
  #   organism (required): string. Name of organism
  #   sequence_length: string. Intervals as "1:500"
  #   other_search_words: string. collection of search terms in one string. E.g.
  #                       "complete mitochondrial"
  #   gene: string. collection of gene names. E.g. "cox1 cox2"
  #   filename: name of the written fasta file
  #   keep: boolean. if TRUE (default) the sequences are written to R as an object. 
  # Returns: if keep = TRUE: object and file written to directory
  #          if keep = FALSE: file written to directory.
  require(rentrez)
  require(purrr)
  # helper function
  paste_special <- function(string, abbreviation){
    
    string <- str_split(string, pattern = " ") %>% unlist %>% 
      paste0(abbreviation)
    if (length(string) > 1){
      string <- paste(string, collapse = " ")}
   return(string)
  }
  
  # Organism
  org <- paste0(organism,"[ORGN]")
  
  # Sequence Length
  if (!is_empty(sequence_length)) {
    sequence_length <- paste_special(sequence_length, "[SLEN]")}
  
  # Gene
  if (!is_empty(gene)) {
    gene <- paste_special(gene, "[GENE]")}
  
  # WORDS
  if(!is_empty(other_search_words)) {
    other_search_words <- paste_special(other_search_words, "[WORD]")}
  
  query <- paste(org,
                 if (!is_empty(sequence_length)) sequence_length,
                 if(!is_empty(other_search_words)) other_search_words,
                 if(!is_empty(gene)) gene)
  
  query <- str_replace_all(query, "   ", " AND ")
  print(query)
  search <- entrez_search(db=database, term=query, retmax = max_entries)
  print(search)
  
  # Extract ids from search
  ids <- search$ids
  
  # fetch sequences
  sequences_fasta <- entrez_fetch(db = database, id = ids, rettype = datatype)
  write(sequences_fasta, 
        file = if (is_empty(filename)){
          paste(paste0(organism, collapse = "_"), datatype, sep = ".")
        } else {
          filename})
}

# Examples
# extract_sequences_genbank("Apis", sequence_length = "1:10000", gene = "COX1", filename = "Apis_cox1")
# 
# extract_sequences_genbank("Ranunculaceae", sequence_length = "100000:200000", other_search_words = "complete", filename = "Ranunculaceae_plastids")
