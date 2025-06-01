library(gifski)

setwd("D:/ngs storage/Natalie/R255X E488QD/emergeR/")
here::i_am("eda/edits/tcc_motifs/tcc_6to8/tcc_6to8_slidingkmer_logos.R")

outdir  <- "eda/edits/tcc_motifs/tcc_6to8/"
outfile <- "tcc_6to8.gif"


# ============================================================================ #
# Stiches together images in the same directory into an animation

frames <- list.files(here(outdir),
                     pattern="\\.png$", full.name=TRUE)
frames <- sort(frames)

n_hold_frames <- 1 # arbitrary
frames <- c(rep(head(frames,1), n_hold_frames),
            frames,
            rep(tail(frames,1), n_hold_frames)
)

gifski(png_files = frames,
       gif_file = here(outdir, outfile),
       width=600, height=600, delay=0.3)
