#' Clean an EMERGe dataset after HT Stream preprocessing.
#' 
#' Reads a FASTQ file and returns a standardized tibble.
#' 
#' @param fastq_infile (String) Path to a FASTQ file.
#' @param template Character(1) DNA template sequence.
#' @param flanks Character(4) flank subsequences in expected order,
#'               (f1, f2, f3, f4).
#' @return A tibble with N10-level summaries: number of reads "n," 
#'         number of edits "k," and other codons "other."
#' @export
#' 
#' @examples 
#' \dontrun{
#' # R255XE488Q. Target is GAA edited to GGA. We take 6-mer flanks around GAA 
#' # and the N10. The full RNA hairpin follows:
#' # GGGAATGGCGCTTCCCGGCA GGAAGT GAA AAGCTG AGGCCGACTTTCTTTCTTTCGTCG GCCTCA 
#' # NNNNNNNNNN TCCTGC CGGGTTGTAGGGTC
#' # Flanks, demarcated above by spaces: GGAAGT AAGCTG GCCTCA TCCTGC
#' template = paste0(
#'   "GGGAATGGCGCTTCCCGGCAGGAAGTGAAAAGCTGAGGCCGACTTTCTTTCTTTCGTCGGCCTCA",
#'   "NNNNNNNNNNTCCTGCCGGGTTGTAGGGTC"
#' )
#' flanks = c("GGAAGT", "AAGCTG", "GCCTCA", "TCCTGC")
#' emergeR::clean(
#'   fastq_infile = "example.fastq", template = template, flanks = flanks
#' )
#' }
#' 
#' @importFrom ShortRead readFastq
#' @importFrom Biostrings DNAString DNAStringSet matchPattern sread extractAt 
#' @importFrom Biostrings subseq xscat
#' @importFrom BiocGenerics width
#' @importFrom IRanges IRanges IRangesList start end 
#' @importFrom tibble tibble
#' @importFrom dplyr count group_by summarise ungroup mutate
#' @importFrom magrittr %>%

clean <- function(
    fastq_infile, 
    template,
    flanks # c(f1, f2, f3, f4)
) {
  size_cdna  <- nchar(template)
  size_left  <- size_cdna - 2L
  size_right <- size_cdna + 1L
  
  reads <- readFastq(fastq_infile)
  reads <- sread(reads[width(reads) %in% c(size_left:size_right)])
  
  flank_set <- DNAStringSet(flanks)
  regions   <- derive_regions(tpl_chr = template, flanks = flanks)
  pieces    <- extractAt(reads, IRangesList(regions$flanks_range))
  
  flat      <- unlist(pieces, use.names = FALSE)
  cmp       <- flat == rep(flank_set, length(reads))
  mat       <- matrix(cmp, nrow = length(flank_set))
  perfect   <- colSums(mat) == length(flank_set)
  
  key_seqs <- as.character(
    subseq(reads[perfect], start(regions$edit_range), end(regions$edit_range))
  )
  n10_seqs <- as.character(
    subseq(reads[perfect], start(regions$n10_range), end(regions$n10_range))
  )
  
  key_seqs <- toupper(key_seqs)
  stopifnot(all(nchar(key_seqs) == 3L))
  
  out <- tibble(n10 = n10_seqs, codon = key_seqs) %>%
    mutate(mid = substr(codon, 2, 2)) %>%
    count(n10, mid, name = "count") %>%
    group_by(n10) %>%
    summarise(
      n = sum(count[mid %in% c("G", "A")]),
      k = sum(count[mid == "G"]),
      .groups = "drop"
    ) %>%
    ungroup()
  
  print(out)
  
  return(out)
}

#' Derives the indices for provided flanks within a DNA template.
#' 
#' @param tpl_chr (string) DNA template used in EMERGe screen, with 'NNNNNNNNNN'
#'                filling the N10 region.
#' @param flanks (string[4]) flanks used to identify successful reads.
#'               See clean() for more information.
#' @return IRanges indices for flanks within the template DNA.
#' 
#' @keywords internal
#' @noRd
#'
#' @importFrom Biostrings DNAString DNAStringSet matchPattern
#' @importFrom IRanges IRanges start end
#' @importFrom S4Vectors xscat
derive_regions <- function(tpl_chr, flanks) {
  stopifnot(length(flanks) == 4L)
  tpl <- DNAString(tpl_chr)
  
  hits <- lapply(flanks, function(f) {
    matchPattern(DNAString(f), tpl, max.mismatch = 0)
  })
  if (!all(lengths(hits) == 1L))
    stop("Flanks not uniquely found in template.")
  
  starts <- vapply(hits, start, 0L)
  ends   <- vapply(hits, end, 0L)
  
  ord    <- order(starts)
  starts <- starts[ord]; ends <- ends[ord]; flanks <- flanks[ord]
  
  flanks_range <- IRanges(start = starts, end = ends)
  
  edit_range  <- IRanges(start = end(flanks_range[1]) + 1L,
                        end   = start(flanks_range[2]) - 1L)
  
  n10_range   <- IRanges(start = end(flanks_range[3]) + 1L,
                        end   = start(flanks_range[4]) - 1L)
  
  list(
    tpl          = tpl,
    flanks_range = flanks_range,
    edit_range   = edit_range,
    n10_range    = n10_range,
    expected     = xscat(DNAStringSet(flanks))
  )
}