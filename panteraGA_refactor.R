#!/usr/bin/env Rscript

# First version using FastGA data.
# Improved stats and recovery of TEs
# Version starting directly from an alignment of two haplotypes
# Cleaned code. Moved to 0.5. Removed surplus options.
# Special version for the VGP project.
# Updates to filter 
# 0.5.2 fix for segment dups
# 0.6.0 Name changed.
# 0.6.1 Included new method to find the edges. TSDs returned. Consensus includes
# now N in positions with no consensus but enough saturation.

pantera_version <- "0.6.1"
options(warn = 0)


###
### Find the path of the script and create an Rlibs folder to store the 
### required libraries if they are missing.
###

stub <- function() {}
thisPath <- function() {  # Path of the script before this.path (Add source!)
cmdArgs <- commandArgs(trailingOnly = FALSE)
if (length(grep("^--file=", cmdArgs)) > -1) {
    # Rscript/R console option
    scriptPath <- normalizePath(dirname(sub("^--file=", "",
                                            cmdArgs[grep("^--file=",
                                                         cmdArgs)])))[1]
  } else if (length(grep("^-f$", cmdArgs)) > 0) {
    # R console option
    normalizePath(dirname(cmdArgs[grep("^-f", cmdArgs) + 1]))[1]
  }  else if (Sys.getenv("RSTUDIO") == "1") {
    # RStudio
    dirname(rstudioapi::getSourceEditorContext()$path)
  } else if (is.null(attr(stub, "srcref")) == FALSE) {
    # 'source'd via R console
    dirname(normalizePath(attr(attr(stub, "srcref"), "srcfile")$filename))
  } else {
    stop("Cannot find file path")
  }
}

### Uncomment if script cannot write to default R libs folder
# if (!file.exists(file.path(thisPath(), "Rlibs"))) {
#   dir.create(file.path(thisPath(), "Rlibs"))
# }
# local({r <- getOption("repos")
# r["CRAN"] <- "https://cran.r-project.org"
# options(repos = r)
# })
# .libPaths(c(file.path(thisPath(), "Rlibs"), .libPaths()))



###
### Auxiliary functions
###

# Prints log messages
lx <- function(x) {
  #  g <- gc()
  cat(paste0(
    format(Sys.time(), "%y-%m-%d:%H:%M:%S"),
    " \u001b[95;1;1m[pantera ",
    pantera_version, 
    #    "]\u001b[0m Mem: ",g[1,2]+g[2,2]," Mb : ",
    "]\u001b[0m ",
    x,
    "\n"
  ))
}

# Installs the required libraries. Probably I should include some messages.
get_libs <- function() {

  if (suppressPackageStartupMessages(!require("this.path", quietly = TRUE))) {
    lx("Intalling package [this.path]")
    install.packages("this.path")
  }
  if (suppressPackageStartupMessages(!require("getopt", quietly = TRUE))) {
    lx("Intalling package [getopt]")
    install.packages("getopt")
  }
  if (suppressPackageStartupMessages(!require("parallel", quietly = TRUE))) {
    lx("Intalling package [parallel]")
    install.packages("parallel")
  }
  if (suppressPackageStartupMessages(!require("ips", quietly = TRUE))) {
    lx("Intalling package [ips]")
    install.packages("ips")
  }
  if (suppressPackageStartupMessages(!require("BiocManager", quietly = TRUE))) {
    lx("Intalling package [BiocManager]")
    install.packages("BiocManager")
  }
  if (suppressPackageStartupMessages(!require("Biostrings", quietly = TRUE))) {
    lx("Intalling package [Biostrings]")
    BiocManager::install("Biostrings")
  }
  if (suppressPackageStartupMessages(!require("DECIPHER", quietly = TRUE))) {
    lx("Intalling package [DECIPHER]")
    BiocManager::install("DECIPHER")
  } 
  if (suppressPackageStartupMessages(!require("bioseq", quietly = TRUE))) {
    lx("Intalling package [bioseq]")
    BiocManager::install("bioseq")
  }
  if (suppressPackageStartupMessages(!require("ape", quietly = TRUE))) {
    lx("Intalling package [ape]")
    install.packages("ape")
  }
  if (suppressPackageStartupMessages(!require("LncFinder", quietly = TRUE))) {
    lx("Intalling package [LncFinder]")
    install.packages("LncFinder")
  }
  if (suppressPackageStartupMessages(!require("seqinr", quietly = TRUE))) {
    lx("Intalling package [seqinr]")
    install.packages("seqinr")
  }
  if (suppressPackageStartupMessages(!require("purrr", quietly = TRUE))) {
    lx("Intalling package [purrr]")
    install.packages("purrr")
  }
  if (suppressPackageStartupMessages(!require("dplyr", quietly = TRUE))) {
    lx("Intalling package [dplyr]")
    install.packages("dplyr")
  }
  if (suppressPackageStartupMessages(!require("data.table", quietly = TRUE))) {
    lx("Intalling package [data.table]")
    install.packages("data.table")
  }
  if (suppressPackageStartupMessages(!require("xgboost", quietly = TRUE))) {
    lx("Intalling package [xgboost]")
    install.packages("xgboost")
  }
  if (!require("stringr", quietly = TRUE)) {
    lx("Intalling package [stringr]")
    install.packages("stringr")
  }
  if (!require("stringi", quietly = TRUE)) {
    lx("Intalling package [stringi]")
    install.packages("stringi")
  }
 
}

# Reads parameters
read_pars <- function() {
  
  spec <- matrix(
    c("sv_file",        "1", 1, "character",   # 1aln input file [required]
      "output_folder",  "o", 1, "character", # Output foder [pantera_output]
      "lib_name",       "b", 1, "character", # Name to append to the output sequences [pantera]
      "threads",        "T", 1, "integer",   # Number of threads [all]
      "min_size",       "s", 1, "integer",   # Min size of polymorphic sequences to investigate [200]
      "max_size",       "l", 1, "integer",   # Max size of polymorphic sequences to investigate [20000]
      "identity",       "i", 1, "double",    # Cutoff (as distance) to cluster in first pass [0.03]
      "identity2",      "y", 1, "double",    # Cutoff (as distance) to cluster in second pass [0.05]
      "min_cl",         "m", 1, "integer",   # Min number of sequences required to cluster [3]
      "Ns",             "n", 1, "integer",   # Max % of Ns allowed in a segment [0]
      "pA_bases",       "p", 1, "integer",   # Min number of bases of a polyA [10]
      "cl_size",        "u", 1, "integer",   # Max number of sequences to cluster [200] 
      "flanking",       "f", 1, "integer",   # Min length of flanking sequences [4000]
      "verbose",        "v", 0, "logical",   # Show log messages
      "log_file",       "r", 0, "logical",   # Create log file 
      "debug",          "d", 0, "logical",   # Keep intermediate files.
      "help",           "h", 0, "logical"),  # Show this help 
    byrow = TRUE, ncol = 4)
  
  opt <- getopt(spec)
  if (!is.null(opt$help)) {
    cat(gsub("\\] \\[","\\]\n\\[",getopt(spec, usage = TRUE)))
    q(status = 1)
  }
  
  if (is.null(opt$log_file)) {
    opt$log_file <- FALSE
  }
  
  if (is.null(opt$verbose)) {
    opt$verbose <- FALSE
  }
  
  if (is.null(opt$debug)) {
    opt$debug <- FALSE
  }
    if (is.null(opt$sv_file)) {
    stop("[--svfile|-1] missing")
   }
  
  if (is.null(opt$output_folder)) {
    opt$output_folder <- "pantera_output"
  }
  
  if (is.null(opt$lib_name)) {
    opt$lib_name <- "pantera"
  }
  
  if (is.null(opt$threads)) {
    opt$threads <- detectCores()
  } else {
    opt$threads <- min(opt$threads, detectCores())
  }
  
  if (is.null(opt$min_size)) {
    opt$min_size <- 100
  }
  
  if (is.null(opt$max_size)) {
    opt$max_size <- 30000
  }
  
  if (is.null(opt$identity)) {
    opt$identity <- 0.03
  }
  
  if (is.null(opt$identity2)) {
    opt$identity2 <- 0.05
  }
  
  if (is.null(opt$min_cl)) {
    opt$min_cl <- 3
  }
  
  if (is.null(opt$Ns)) {
    # Maximum percentage of Ns in segment.
    opt$Ns <- 0
  }
  
  if (is.null(opt$pas)) {
    # Maximum percentage of Ns in segment.
    opt$pas <- 10
  }
  
  if (is.null(opt$pag)) {
    # Maximum percentage of Ns in segment.
    opt$pag <- 5
  }
  
  if (is.null(opt$cl_size)) {
    # Maximum size of cluster to process
    opt$cl_size <- 300
  }
  
  if (is.null(opt$flanking)) {
    # min length of mapping flanking sequences to define an SV
    opt$flanking <- 1000
  }
  return(opt)
}

# Reads a fasta file
ffasta <- function(f) {
  check <- file.size(f)
  if (check == 0 | is.na(check)) {
    return(data.table(name = as.character(), seq = as.character()))
  } else {
  fa_raw <- fread(f, header = F, fill = T, sep = "\n")
  fa_raw[, h := grepl(">", V1)]
  fa_oneline <- fa_raw[, .(paste0(V1, collapse = "")), by = rleid(h)]
  return(data.table(name = fa_oneline[rep(c(TRUE, FALSE), length = .N)]$V1, 
                    seq = fa_oneline[rep(c(FALSE, TRUE), length = .N)]$V1))
  }
}

# Reads the two ONEcode files extracted from a 1aln file with svfind
ffastasv <- function(f) {
  check <- file.size(f)
  if (check == 0 | is.na(check)) {
    return(data.table(name = as.character(), seq = as.character()))
  } else {
    fa_raw <- fread(cmd=paste("seqconvert -fa",f, " 2> /dev/null"),
                    header = F,
                    fill = T,
                    sep = "\n")
    fa_raw[, h := grepl(">", V1)]
    fa_oneline <- fa_raw[, .(paste0(V1, collapse = "")), by = rleid(h)]
    return(data.table(name = fa_oneline[rep(c(TRUE, FALSE), length = .N)]$V1, 
                      seq = fa_oneline[rep(c(FALSE, TRUE), length = .N)]$V1))
  }
}

parseONEview <- function(f) {
  lx(paste("Procesing file =", f))
  check <- file.size(f)
  if (check == 0 | is.na(check)) {
    end_pantera(".1aln file missing")
  } else {
  system(paste0("svfind -x 12 -f ", opt$flanking, " -a ",opt$lib_name,"hapa -b ",opt$lib_name,"hapb ", f, " > /dev/null 2>&1"))
    ova <- fread(cmd = paste0("ONEview ",opt$lib_name,"hapa"), fill = T)
    ovb <- fread(cmd = paste0("ONEview ",opt$lib_name,"hapb"), fill = T)
    ova[,f:="a"]
    ovb[,f:="b"]
    ov <- rbind(ova,ovb)
    sequences <- ov[V1=="S",2:3][,V2:=as.numeric(V2)]
    file <- ov[V1=="S",f]
    overhangs <- ov[V1=="Q",2:3][,V2:=as.numeric(V2)]
    flanks <- ov[V1=="F",2:3][,V2:=as.numeric(V2)][,V3:=as.numeric(V3)]
    name <- ov[V1=="I",3]
    variant <- ov[V1=="V",2:4][,V2:=as.numeric(V2)][,V3:=as.numeric(V3)][,V4:=as.numeric(V4)]
    source <- ov[V1=="B",2:4][,V2:=as.numeric(V2)][,V3:=as.numeric(V3)][,V4:=as.numeric(V4)]
    ovall <- cbind(file,name,variant[,1:3], source[,1:3],flanks[,1:2], sequences[,1:2], overhangs[,1:2])
    colnames(ovall) <- c("file","name","v_name", "v_start","v_end","s_name", "s_start","s_end", "flank_l", "flank_r", "seq_len","seq","ov_len","overhang")
    lx(paste("Total segments in 1aln =", nrow(ovall)))
    ## Remove segments that map to the same spot, we remove the last digit to account for small missmatches in the breakpoint
    ovall[,ins_point :=paste0(file,"-",s_name,":",substr(s_start,1,nchar(s_start)-1),"-",substr(s_end,1,nchar(s_end)-1))]
    ovall <- ovall[!(ins_point %in% ovall$ins_point[duplicated(ovall$ins_point)]) ]  
    ovall <- ovall[seq_len >= opt$min_size & seq_len <= opt$max_size]
    lx(paste("Total segments filtered =", nrow(ovall)))
  ### It seems that removing duplicated insertions is enough to improve the data, the flanking segments distribution improves a lot, for now we do not do anything else with them  
  #  flankcut <-  quantile(c(ovall$flank_r,ovall$flank_l))
  #  ovcut <- ovall[flank_l>flankcut[2] & flank_r>flankcut[2]]
    ### Remove all but one copy when the same segment can map to several places.
    ovall <-ovall[order(file,v_name,v_start)]
    ovall[,ov:=shift(v_start,-1), by=v_name]
    ovall[,ov:=ov-v_end]
    ovall[is.na(ov) ,ov:=1]
    ovall[ov < -seq_len ,ov:=1]
    ovall <- ovall[ov>0]
   
    unlink(paste0(opt$lib_name,"hapa"))
    unlink(paste0(opt$lib_name,"hapb"))
    
    if (opt$Ns > 0) {
      ovall[, Ns := nchar(ovall$seq) - nchar(gsub("N", "", ovall$seq))]
      ovall <- ovall[Ns <= (len * opt$Ns)]
      lx(paste("Number of valid Ns segments =", nrow(ovall)))
    }
    ovall[,name:=paste0(">",.I)]
  }
  return(ovall[,c("name","seq")])
}

# Reverse complement of a DNA sequence
rc <- function(s) {
  stri_reverse(chartr("ACGTacgt","TGCAtgca",s))
}

# Writes a fasta file
wfasta <- function(fasta_data, f) {
  fileconn <- file(f)
  writeLines(t(fasta_data), fileconn)
  close(fileconn)
}

# End and show message
end_pantera <- function(message) {
  cat(paste0(
    format(Sys.time(), "%y-%m-%d:%H:%M:%S"),
    " \u001b[95;1;1m[pantera ",
    pantera_version,
    "]\u001b[0m ",
    "\u001b[33;101;1m",
    "\u2639 ",
    "\u001b[0m : ",
    message,
    "\n"
  ))
  quit(save = "no")
}

extract_all_substrings <- function(string) {
  n <- stri_length(string)
  
  # Generate all combinations of start positions and lengths
  combinations <- do.call(rbind, lapply(1:n, function(len) {
    max_start <- n - len + 1
    if (max_start >= 1) {
      data.frame(start = 1:max_start, end = 1:max_start + len - 1)
    } else {
      data.frame(start = integer(0), end = integer(0))
    }
  }))
  
  # Extract all substrings using vectorized stri_sub
  stri_sub(string, combinations$start, combinations$end)
}

longest_common_substrings <- function(string1, string2) {
  subs1 <- extract_all_substrings(string1)
  subs2 <- extract_all_substrings(string2)
  
  # Find common substrings
  common_subs <- intersect(subs1, subs2)
  
  if (length(common_subs) == 0) {
    return(character(0))
  }
  
  # Find maximum length
  max_length <- max(stri_length(common_subs))
  
  # Return all substrings with maximum length
  common_subs[stri_length(common_subs) == max_length]
}

# Converts from dna to binDNA format
reformatDNA <- function(dna) {
  temp <- matrix(as.character(dna),
                 nrow = (length(row.names(dna))),
                 dimnames = dimnames(dna))
  temp <- apply(temp, 1, function(x) {
    paste0(x, collapse = "")
  })
  return(temp)
}


###
### Main functions 
###

# Init function
init_pantera <- function() {
  
  dir.create(opt$output_folder,
             showWarnings = FALSE,
             recursive = T)
  if (opt$log_file) {
  sink(paste0(opt$output_folder, "/pantera.log"), split = opt$verbose) 
    }
lx(paste("pantera", pantera_version, "\n"))

## Confirm dependencies are available
if (!(nchar(Sys.which("mafft")) > 0)) {
  end_pantera("mafft not found")
}
if (!(nchar(Sys.which("getorf")) > 0)) {
  end_pantera("getorf not found")
}
if (!(nchar(Sys.which("blastn")) > 0)) {
  end_pantera("blastn not found")
}

if (!file.exists(file.path(this.dir(), "model/xgbmodel_typesnames_fix"))) {
  end_pantera("Model not found.")
}

lx(paste("pantera path:", this.path()))
lx(paste("Gfas list:", opt$gfas_list))
lx(paste("Output:", opt$output_folder))
lx(paste("Threads to use:", opt$threads))
lx(paste("Cores detected:", detectCores()))
lx(paste("Min. size:", opt$min_size))
lx(paste("Max. size:", opt$max_size))
lx(paste("Identity for clustering:", opt$identity))
lx(paste("Identity for clustering 2:", opt$identity2))
lx(paste("Min. sequences to create a consensus:", opt$min_cl))
lx(paste("Max percentage of Ns:", opt$Ns))
lx(paste("Max. size of cluster:", opt$cl_size))
lx(paste("Minimum size of flanking sequences:", opt$flanking))

return(0)
}

# Reads polymorfic segments from an alignment
get_segments <- function(segments_unique) {
  ## Extract unique segments
  lx(paste("Procesing file =", opt$sv_file))
  system(paste0("svfind -f ", opt$flanking, " -a ",opt$lib_name,"hapa -b ",opt$lib_name,"hapb ", opt$sv_file, " > /dev/null 2>&1"))
  segmentsa <- ffastasv(paste0(opt$lib_name,"hapa"))
  segmentsb <- ffastasv(paste0(opt$lib_name,"hapb"))
  unlink(paste0(opt$lib_name,"hapa"))
  unlink(paste0(opt$lib_name,"hapb"))
  segments <- rbind(segmentsa,segmentsb)
  segments[, len := nchar(seq)]
  segments <- segments[len >= opt$min_size & len <= opt$max_size]
  segments[,b:=gsub(".*_","",name)]
  segments <- segments[order(-len)][,.SD[1],b]  #[,.(name,seq)]
  segments[,a:=gsub("_.*","",name)]
  segments <- segments[,.SD[1],a]# [,.(name,seq)]
  segments[,ca:=gsub(">","",gsub(":.*","",a))]
  segments[,sa:=as.numeric(gsub(".*:","",gsub("-.*","",a)))]
  segments[,ea:=as.numeric(gsub(".*-","",a))]
  segments[,cb:=gsub(":.*","",b)]
  segments[,sb:=as.numeric(gsub(".*:","",gsub("-.*","",b)))]
  segments[,eb:=as.numeric(gsub(".*-","",b))]
  segments <-segments[order(ca,sa)]
  segments[,ov:=shift(ea,-1), by=ca]
  segments[,ov:=ov-ea]
  segments[,rn:=.I]
  dups <-segments[segments$ov<0 | segments[,shift(ov,1), by=ca]$V1<0]
  ndups <- nrow(dups)
  dups$g <-  cumsum(c(1L,diff(dups$rn)>1))
  segments <- segments[!(name %in% dups$name)]
  dups <- dups[dups[, .I[len == max(len)], by=g]$V1]
  lx(paste("Removed overlaps =", ndups - nrow(dups)))
  segments <- rbind(segments,dups[,1:13])
  segments <- segments[,name:=paste0(">",basename(opt$sv_file),
                                     "_",gsub(">","",a), collapse = ""), by =.I]
  lx(paste("Number of valid size segments =", nrow(segments)))
  if (opt$Ns > 0) {
    segments[, Ns := nchar(segments$seq) - nchar(gsub("N", "", segments$seq))]
    segments <- segments[Ns <= (len * opt$Ns)]
    lx(paste("Number of valid Ns segments =", nrow(segments)))
    segments$Ns <- NULL
  }
  segments_unique <- segments[,.(name,seq)]
  lx(paste("TOTAL number of segments =", nrow(segments_unique)))
  return(segments_unique)
}

# Obtains and stores the polymorphic segments
read_poly <- function() {
  if (!file.exists(paste0("all_segments.fa"))) {
    segments_unique <- parseONEview(opt$sv_file)
    segments_unique[,seq:=toupper(seq)]
    wfasta(segments_unique[, c("name", "seq")],
           paste0(opt$output_folder, "/all_segments.fa"))
  }
  segments_unique[,len:=nchar(seq)]
return(segments_unique)  
}  

# First pass, finds repeats but does not try to generate any consensus.  
find_repeats <- function() {
  process_zone <- function(zone) {
    start <- as.numeric(zone[1]) - zones_interval_overlap
    end <- as.numeric(zone[2]) + zones_interval_overlap
    discards <- 0
    if (start != end) {
      # If they have the same value they were processed in another zone
      if (nrow(segments_unique[len >= start][len <= end])> 0) {
        segment_u <- segments_unique[len >= start][len <= end]
        segment_u <- segment_u[!duplicated(segment_u)]
        if (nrow(segment_u) >= opt$min_cl)  {
          lx(paste0("Num. of segments zone ",start," - ",end," : ", 
                    nrow(segment_u)))
          segment_sets <- split(segment_u, 
                                rep(1:(nrow(segment_u) %/% opt$cl_size +1), 
                                length.out = nrow(segment_u), 
                                each = ceiling(nrow(segment_u) / 
                                               (nrow(segment_u) %/% 
                                                opt$cl_size +1))))
          candidates <- data.table(name = as.character(), 
                                   seq = as.character())
          for (ss in 1:length(segment_sets)) {
            lx(paste0("Num. of segments zone ", start," - ", 
                      end," : cluster : ",ss))
            sg <- segment_sets[[ss]]
            sgrc <- data.table(name=sg$name, seq=rc(sg$seq))
            sgrc[,name:=paste0(name,"_RC")]
            sg_seq <- rbindlist(list(sg[,c("name","seq")],sgrc))
            set2clus <- DNAStringSet(sg_seq$seq)
            list_clust <- data.table(Clusterize(set2clus, 
                                                cutoff = opt$identity, 
                                                verbose = F, 
                                                penalizeGapLetterMatches = T, 
                                                includeTerminalGaps = T,
                                                processors = 1, 
                                                minCoverage = 0.9))
            list_clust[,name:= sg_seq$name]
            good_clust <- list_clust[,.N,cluster][N>=opt$min_cl]
          
            lx(paste("# of clusters ", start, "-", end, ":", 
                     length(unique(good_clust))))
         
            candidates <- rbindlist(list(candidates, 
                  sg[name %in% list_clust[cluster %in% good_clust$cluster]$name, 
                     c("name", "seq")]))
          }
          if (nrow(candidates) > 0) {
            wfasta(candidates, paste0("candidates_", start, "_", end, ".fa"))
          }
        }
      }
    }
    return(0)
  }
  
  round <- 1
  lx(paste("Reading data for round:", round))
  segments_unique <- ffasta("all_segments.fa")
  lx(paste("TOTAL Unique segments:", nrow(segments_unique)))
  segments_unique <- segments_unique[!duplicated(segments_unique)]
  segments_unique[, len := nchar(seq)]
  segments_unique <- segments_unique[order(-len)]
  lx(paste("TOTAL Unique segments no dups:", nrow(segments_unique)))
  lx(paste("Starting loop:", round))
  lx(paste("Processing:", nrow(segments_unique), "segments"))
  lx(paste("Largest segment:", max(segments_unique$len)))
  lx(paste("Smallest segment:", min(segments_unique$len)))
  dir.create(paste0("loop_", round), showWarnings = FALSE)
  setwd(paste0("loop_", round))
  z <-
    seq(
      min(segments_unique$len),
      max(segments_unique$len) + zones_interval,
      zones_interval
    )
  if (length(z) < 2) {
    zones <-
      list(data.frame(
        start = min(segments_unique$len),
        end = max(segments_unique$len)
      ))
  } else {
    zones <- asplit(data.frame(start = z[-length(z)], end = z[-1]), 1)
  }
  
  lx(paste("Processing:", length(zones), "windows"))
  loop_exit <- mclapply(
    rev(zones),
    process_zone,
    mc.preschedule = FALSE,
    mc.cores = opt$threads
  )

  lx(paste("Total sequences discarded:", sum(as.numeric(unlist(loop_exit)))))
  system(paste0("mv ../all_segments.fa ../pantera_lib_", round - 1, ".fa"))
  system("cat candidates*.fa > ../segments_candidates.fa")
  setwd("..")
}

# Second pass. Cluster and generate consensus. 
cluster_results <- function() {

  process_zone_2 <- function(zone) {
    start <- as.numeric(zone[1]) - zones_interval_overlap
    end <- as.numeric(zone[2]) + zones_interval_overlap
    discards <- 0
    if (start != end) {
      lx(paste("Procesing segments:", start, "-", end))
      if (nrow(segments_unique[len >= start][len <= end])> 0) {
        segment_u <- segments_unique[len >= start][len <= end]
        segment_u <- segment_u[!duplicated(segment_u)]
        lx(paste("Number of sequences:", nrow(segment_u)))
          if (nrow(segment_u) >= opt$min_cl)  {
            segment_sets <- split(segment_u, 
                                  rep(1:(nrow(segment_u) %/% opt$cl_size +1), 
                                      length.out = nrow(segment_u), 
                                      each = ceiling(nrow(segment_u) / 
                                                     (nrow(segment_u) %/%
                                                      opt$cl_size +1))))
            lx(paste("Segment sets", start, "-", end, ":", 
                     length(segment_sets)))
            for (ss in 1:length(segment_sets)) {
              sg <- segment_sets[[ss]]
              sgrc <- data.table(name=sg$name, seq=rc(sg$seq))
              sgrc[,name:=paste0(name,"_RC")]
              sg_seq <- rbindlist(list(sg[,c("name","seq")],sgrc))
              set2clus <- DNAStringSet(sg_seq$seq)
              list_clust <- data.table(Clusterize(set2clus, 
                                                  cutoff = opt$identity2, 
                                                  verbose = F, 
                                                  penalizeGapLetterMatches = T, 
                                                  includeTerminalGaps = T,
                                                  processors = 1, 
                                                  minCoverage = -0.9))
              lx(paste("Segment sets", start, "-", end, ": clusters", 
                       length(unique(list_clust$cluster))))
              list_clust[,name:= sg_seq$name]
              list_clust[,name:=gsub("_RC","",name)]
              clusters <- list_clust[order(name)][,.(list(name)), cluster]
              clusters <- clusters[!duplicated(clusters$V1)]
              clusters[,n:=length(V1[[1]]), by=.I]
              clusters <- clusters[clusters$n>=opt$min_cl,] 
              consensi <- data.table(name = as.character(), 
                                     seq = as.character())
          if ( nrow(clusters) > 0) {
            for (u in unique(clusters$cluster)) {
              seqs_clust <- sg[name %in% unlist(clusters[cluster==u]$V1),]$seq
              clust_temp <- strsplit(seqs_clust, "")
              names(clust_temp) <- sg[name %in% 
                                      unlist(clusters[cluster==u]$V1),]$name
              seqs <- as.DNAbin(clust_temp)
              if (length(seqs) > 1) {
                ali <- ips::mafft(seqs,
                             #      method = "globalpair",
                             #      maxiterate = 2,
                             options = c("--adjustdirection"),
                             ep = 0.123,
                             thread = -1,
                             exec = "mafft"
                )
                alim <- reformatDNA(ali)
                cons <- toupper(consensusString(DNAStringSet(alim), 
                                                ambiguityMap = "-", 
                                                threshold = cons_threshold))
                
                ### Exploring how to find the TSDs. The cons gives as already the positions. Maybe get pos +-1,checkt that one, and then extend to both sides?
                
                cMatrix <- consensusMatrix(DNAStringSet(alim))
                dtMatrix <- as.data.table(matrix(unlist(cMatrix), nrow = nrow(cMatrix)))
                saturation <- unlist(lapply(dtMatrix[1:4],sum))/length(seqs)
                conserv <- unlist(lapply(dtMatrix[1:4],max))/length(seqs)
               # kclus <- kmeans(conserv,2)
                consstring <- paste0(as.numeric(conserv>=cons_threshold+.4),collapse = "")
               # good <- which(kclus$size == max(kclus$size), kclus$size)
               # consstring <- paste0(kclus$cluster,collapse ="")
                cons_s <- stri_locate_first_fixed(consstring,paste0(rep(1,max(2,(12-floor(sqrt(nrow(ali)-2))))),collapse =""))[1]
                cons_e <- stri_locate_last_fixed(consstring,paste0(rep(1,max(2,(12-floor(sqrt(nrow(ali)-2))))),collapse =""))[2]
                edges <- data.table(l=stri_reverse(gsub("-","",substr(alim,1,cons_s-1))), r=gsub("-","",substr(alim,cons_e+1,nchar(alim[1]))))
                tsds <- unlist(lapply(2:13,function(x){sum(diag(adist(substr(edges$r,1,x),stri_reverse(substr(edges$l,1,x)))))}))
                tsd_len <- min(which(tsds==min(tsds),tsds))+1
                tsds_conf <- round(1- tsds[tsd_len-1]/(2*length(tsds)),2)
                tsds_motif <- consensusString(DNAStringSet(c(substr(edges$r,1,tsd_len),stri_reverse(substr(edges$l,1,tsd_len)))), ambiguityMap=IUPAC_CODE_MAP,
                                threshold=0.25, shift=0L, width=NULL)
                ### End of exploration
              
                
              } else {
                cons <- toupper(paste0(unlist(as.character(seqs)), 
                                       collapse = ""))
              }
              cons <- gsub("-","N",cons)
              cons <- substr(cons,cons_s,cons_e)
              cons <- paste0(strsplit(cons,"")[[1]][saturation[cons_s:cons_e]>saturation_threshold],collapse="")
      
              nam <- paste0(">CONS-", start, "-", end, "-",u, "-", 
                            nchar(cons), "_clus",  length(clust_temp), "_tsdl",tsd_len , "_tsdc",tsds_conf, "_tsdm",tsds_motif,"@@")
              consensi <- rbindlist(list(consensi, data.table(name = nam, 
                                                              seq = cons)))
            }
          } 
        }
                  if (nrow(consensi) > 0) {
                    wfasta(consensi, paste0("consensi_", start, "_",
                                             end, ".fa"))
                  }
      }
      }
    }
    lx(paste("Zone:", start, "completed."))
    return(discards)
  }
  opt$Ns <- 0
  # zones_interval <- 3000
  zones_interval_overlap <- 0
  cons_threshold <- 0.4
  saturation_threshold <- 0.6
  opt$cl_size <- 600
    segments_unique <- ffasta("segments_candidates.fa")
    if (nrow(segments_unique) > 0) {
      lx(paste("TOTAL Unique segments:", nrow(segments_unique)))
      segments_unique <- segments_unique[!duplicated(segments_unique)]
      segments_unique[, len := nchar(seq)]
      segments_unique <- segments_unique[order(-len)]
      lx(paste("TOTAL Unique segments no dups:", nrow(segments_unique)))
      lx(paste("Starting loop 2"))
      lx(paste("Processing:", nrow(segments_unique), "segments"))
      lx(paste("Largest segment:", max(segments_unique$len)))
      lx(paste("Smallest segment:", min(segments_unique$len)))
      dir.create("loop_2", showWarnings = FALSE)
      setwd("loop_2")
      segments_unique[,csum:=cumsum(len)]
      segments_unique[,g:=csum %/% (opt$cl_size * 1000)]
      zones <- segments_unique[,.(min(len),max(len)),g][,2:3]
      zones <- asplit(zones,1)
      # z <- seq(1, nrow(segments_unique), ceiling(nrow(segments_unique) / ceiling(nrow(segments_unique) / opt$cl_size)))
      # if (z[length(z)] != nrow(segments_unique)) {
      #   z <- c(z,nrow(segments_unique))
      # }
      # if (length(z) < 3) {
      #   zones <-
      #     list(data.frame(
      #       start = min(segments_unique$len),
      #       end = max(segments_unique$len)
      #     ))
      # } else {
      #   s = segments_unique[z[-length(z)]]$len +1
      #   e = segments_unique[z[-1]]$len
      #   s[1] <- s[1] - 1
      #   zones <- asplit(data.frame(start = s, end = e), 1)
      # }

      lx(paste("Processing:", length(zones), "windows"))
      loop_exit <- mclapply(
        rev(zones),
        process_zone_2,
        mc.preschedule = FALSE,
        mc.cores = opt$threads
      )
      lx(paste("Total sequences discarded:", 
               sum(as.numeric(unlist(loop_exit)))))
      system("cat consen*.fa >> ../all_consensi.fa")
      setwd("..")
    }
}

# Classify TE models
classify_tes <- function() {
  final <- ffasta("all_consensi.fa")
  final <- final[!duplicated(final$seq)]
  final <- final[order(-nchar(seq))]
  ### Parse information in final to include later
  cluster_n <- as.numeric(gsub("_tsdl.*","",gsub(".*_clus","", final$name)))
  tsd_l <- as.numeric(gsub("_tsdc.*","",gsub(".*_tsdl","", final$name)))
  tsd_c <- as.numeric(gsub("_tsdm.*","",gsub(".*_tsdc","", final$name)))
  tsd_m <- gsub("@@","",gsub(".*_tsdm","", final$name))
  final[, name := paste0(">", opt$lib_name, "_", 1:nrow(final), 
                         "-", gsub(".*-","",name))]
  wfasta(final, paste0(opt$lib_name, "-consensi.fa"))
  lx(paste("Final consensi: ", nrow(final)))
  file <- paste0(opt$lib_name, "-consensi.fa")
  lx("Reading orfs")
  system(
    paste0(
      "getorf -sequence ",
      file,
      " --outseq temp.orfs -minsize 300 &>/dev/null; blastp -num_threads 8 -query temp.orfs -db ",
      this.dir(),
      "/libs/RepeatPeps.lib -outfmt 6 -evalue 1e-10 > orfs.tbl"
    )
  )
  lx("Reading data")
  orfs <- fread("orfs.tbl", header = F)
  if (nrow(orfs) > 0) {
  types <-
    read.table(file.path(this.dir(), "model/xgbmodel_typesnames_fix"),
               sep = "\n")
  types <- types$V1
  features <-
    read.table(file.path(this.dir(), "model/xgbmodel_featurenames"),
               sep = "\n")
  features <- features$V1
  xgb.fit <-
    xgb.load(file.path(this.dir(), "model/xgbmodel_20231030"))
  lx("Processing")
  orfs[, V1 := gsub("_[0-9]*$", "", V1)]
  orfs[, prot := gsub(".*#", "", V2)]
  
  orfs_table <-
    dcast(
      orfs,
      V1 ~ prot,
      value.var = "V12",
      fun.agg = function(x) {
        suppressWarnings(max(x))
      }
    )
  orfs_table[is.na(orfs_table), ] <- 0
  orfs_table[orfs_table == -Inf, ] <- 0
  
  names <- orfs_table$V1
  orfs_table$V1 <- NULL
  
  colnames(orfs_table) <- make.names(colnames(orfs_table))
  if (sum(!(colnames(orfs_table) %in% features)) > 0) {
    orfs_table[, colnames(orfs_table)[!(colnames(orfs_table) %in% features)] :=
                 NULL]
  }
  orfs_table[, features[!(features %in% colnames(orfs_table))] := 0]
  setcolorder(orfs_table, features)
  
  lx("Get predictions")
  xgb.pred2 = predict(xgb.fit, as.matrix(orfs_table), reshape = T)
  xgb.pred2 = as.data.frame(xgb.pred2)
  colnames(xgb.pred2) = types
  xgb.pred2$prediction = apply(xgb.pred2, 1, function(x)
    colnames(xgb.pred2)[which.max(x)])
  xgb.pred2$prob = apply(xgb.pred2[1:ncol(xgb.pred2) - 1], 1, function(x)
    max(x))
  result <-
    data.table(
      Name = names,
      Prediction = xgb.pred2$prediction,
      Probability = xgb.pred2$prob
    )
  fwrite(result, paste0(file, ".predictions"), sep = "\t")
  
  
  result <- fread(paste0(file, ".predictions"))
  result <- result[, name := paste0(">", Name)]
  }
  
  lx("Merge")
  tes <- ffasta(file)
  if (exists("result")) {
    final <- merge(tes, result, all.x = T)
    final[is.na(Prediction), Prediction := "Unknown"]
  } else {
    final <- tes
    final[,Prediction := "Unknown"]
  }
  final[, short_Prediction:=  gsub(".*/","",Prediction)]
  final[,ix:=1:.N, by = short_Prediction]
  final <- final[order(-nchar(seq))]
  final[,tsd_l:=tsd_l]
  final[,tsd_c:=tsd_c]
  final[,tsd_m:=tsd_m]
  final[,cluster:=cluster_n]
  final[, name := paste0(">",short_Prediction,"_",ix,"-",opt$lib_name, 
                         "#", Prediction), by=1:nrow(final)]
  # if (nrow(final[Probability < 0.9]) > 0) {
  # final[Probability < 0.9, name := paste0(name, "_LowConf"), by=.I]
  # }
  final <- final[!duplicated(final$seq)]
 # wfasta(final[, c("name", "seq")], paste0(opt$lib_name, "-pantera-final.fa"))
  return(final[, c("name", "seq", "cluster", "tsd_l","tsd_c","tsd_m")])

}

# Obtain structural stats from TE models and recover or filter some.
stats_tes <- function() {
  lx("Stats start")
  tes <- final[, c("name", "seq", "cluster", "tsd_l","tsd_c","tsd_m")]
  file <- paste0(opt$lib_name, "-consensi.fa")
  temp <- paste0("temp",make.names(gsub(".*/","",file)))
  dir.create(temp)
  setwd(temp)
  lente <- nchar(tes$seq)
  tes$lente <- lente
  three_orfs <- function(seq) {
    orfs <- rbindlist(find_orfs(seq, reverse.strand = TRUE, 
                                max.only = F), fill = T)
    orfs_list <- orfs[,first(ORF.Len), ORF.Stop]
    orfs_list_ordered <- orfs_list[order(-V1)]$V1
    return(data.table(orf1 = orfs_list_ordered[1],
                    orf2 = orfs_list_ordered[2],
                    orf3 = orfs_list_ordered[3]))
  }
  orfs_tes <- rbindlist(mclapply(tes$seq,three_orfs))
  tes <- cbind(tes,orfs_tes)
  wfasta(tes[, c("name","seq")], paste0("tmp-tes"))
  system(paste0("makeblastdb -in tmp-tes -dbtype nucl 1> /dev/null"))
  te_data <- fread(cmd= paste0("blastn -query tmp-tes -db tmp-tes -task blastn -num_threads ", 
                               opt$threads, 
                               " -evalue 5000 -outfmt 6 -word_size 5 -gapopen 4 -gapextend 1 -reward 1 -penalty -1"), header = F)
  if (nrow(te_data) > 0) {
    colnames(te_data) <-c("qseqid","sseqid", "pident" , "length","mismatch", 
                          "gapopen","qstart","qend","sstart","send", 
                          "evalue", "bitscore")
    te_data_mix <- te_data[qseqid != sseqid]
    te_data_rep <- te_data[qseqid == sseqid & length >20 & 
                             (qstart != sstart | qend != send)][,.(N=.N,maxRep=max(length)),qseqid]
    te_data  <- te_data[!(qseqid == sseqid & qstart == sstart & qend == send)]
    te_data <- te_data[qseqid == sseqid]
    te_data[,len:=abs(qend-qstart)]
    te_data <- te_data[order(qstart)][order(-len)]
    te_data <- merge(te_data, 
                     tes[,c("name","lente", "orf1", "orf2", "orf3")]
                     [,name:=gsub(">","",name)], by.x = "qseqid", by.y="name")
    te_data[,lgap:=.(min(qstart,qend,sstart,send)-1), by=1:nrow(te_data)]
    te_data[,rgap:=.(lente-max(qstart,qend,sstart,send)), by=1:nrow(te_data)]
    
    te_data[,check:=((qend+qstart-lente)/2) * ((send+sstart-lente)/2)<0, 
            by=1:nrow(te_data)]  ### Is each match at one side of the center?
    te_data <- te_data[check==T]
    te_data[,check_type:=((qend-qstart)*(send-sstart))<0, by=1:nrow(te_data)]
    te_data[,tgap :=lgap+rgap]
    te_data <- te_data[order(tgap)][,.SD[1],qseqid]
    # te_data_mix <- merge(te_data_mix, 
    #                      tes[,c("name","lente", "orf1", "orf2", "orf3")]
    #                      [,name:=gsub(">","",name)], by.x = "qseqid", 
    #                                                  by.y="name")
    # te_data_mix <- merge(te_data_mix, 
    #                      tes[,c("name","lente", "orf1", "orf2", "orf3")]
    #                      [,name:=gsub(">","",name)], by.x = "sseqid", 
    #                                                  by.y="name")
    # te_data_mix[,lqgap:=.(min(qstart,qend)-1), by=.I]
    # te_data_mix[,rqgap:=.(lente.x-max(qstart,qend)), by=.I]
    # te_data_mix[,lsgap:=.(min(sstart,send)-1), by=.I]
    # te_data_mix[,rsgap:=.(lente.y-max(sstart,send)), by=.I]
    # te_data_mix <- te_data_mix[!grepl("Unknown", sseqid)]
    # te_data_mix <- te_data_mix[lsgap<10 | rsgap<10]
    te_data$type <- "LTR"
    te_data[check_type==T,type:= "TIR"]

   tes <- merge(tes[,name:=gsub(">","",name)], 
                te_data[,c("qseqid","length","lgap","rgap", "type")], 
                         by.x = "name", by.y = "qseqid", all.x=T)
   tes[is.na(type), rgap :=0]
   tes[is.na(type), lgap :=0]
   tes[is.na(type), length :=0]
   
   tes <- merge(tes,te_data_rep, by.x = "name", by.y = "qseqid", all.x=T)
   #  tes[cluster < (opt$min_cl + floor(sqrt(20000 %/% lente)/2) ), pass:=F]
   #  cutoff <- quantile(tes$cluster, 0.8)
   #  tes[cluster > cutoff, pass:= T]
   #  tes[(length > (lente / 2)), pass:=F]
   tes_count <- nrow(tes)
   
   tes$pas <- unlist(lapply(stri_locate_all_regex(tes$seq,paste0(strrep("A",opt$pas),"|",strrep("T",opt$pas))),function(x) {min(unlist(x))-1}))
   tes$eas <- nchar(tes$seq)-unlist(lapply(stri_locate_all_regex(tes$seq,paste0(strrep("A",opt$pas),"|",strrep("T",opt$pas))),function(x) {max(unlist(x))}))
   tes[,pa:=min(pas,eas), by=.I]
   tes[substr(seq,1,5)=="AAAAA" | substr(seq,1,5)=="TTTTT" | substr(seq,lente-4,lente)=="AAAAA" | substr(seq,lente-4,lente)=="TTTTT" ,pa:=0, by = .I]
   
   tes[,pass:=F]
   
   tes[is.na(pa),pa:=1000]
   tes[,TR1:=maxRep>length & maxRep>1000]
   tes[,TR2:=maxRep>lente/2.5]
   tes[is.na(TR1), TR1:=F]
   tes[is.na(TR2), TR2:=F]
   tes[TR1!=T & TR2!=T, pass:=T]
   tes[type == "TIR" & lgap < 3 & rgap< 3, pass:=T]
 
  lx(paste("TRs discards:", nrow(tes[pass==F])))
  }
  setwd("..")
  unlink(temp, recursive = TRUE)
  tes <- tes[order(-lente)]
  tes_fa <- tes[pass==T, c("name", "seq")]
  tes_fa[,name:=paste0(">",name, collapse = ""), .I]
  wfasta(tes_fa, paste0(opt$lib_name, "-pantera-final.fa"))
  tes_dis <- tes[pass==F, c("name", "seq")]
  tes_dis[,name:=paste0(">",name, collapse = ""), .I]
  wfasta(tes_dis, paste0(opt$lib_name, "-pantera-discards.fa"))
  stats_data <- tes[,c("name", "cluster", "tsd_l","tsd_c","tsd_m","pass","lente","type","pa",
                       "length","lgap","rgap","orf1","orf2","orf3","maxRep")]
  colnames(stats_data) <- c("name", "cluster_size", "TSD_length","TSD_confidence","TSD_motif","pass","TE_len", 
                            "struct_type","polyA","struct_len","left_gap",
                            "right_gap","orf1","orf2","orf3","maxTR")
  fwrite(stats_data, paste0(opt$lib_name, "-pantera-final.stats"), 
         quote =  F, row.names = F, sep ="\t")
}

cleanup <- function() {
  
  lx("Library succesfully created.")
  unlink("loop_1", recursive = TRUE)
  unlink("loop_2", recursive = TRUE)
  unlink("orfs.tbl")
  unlink("pantera_lib_0.fa")
  file.rename("segments_candidates.fa", paste0(opt$lib_name, "-segments_candidates.fa"))
  unlink("temp.orfs")
  unlink("*consensi*")
  setwd(current)
  lx(paste0(   "\u001b[32;2m",
               "\u263A",
               "\u001b[0m : ","End of process."))
}


###
### MAIN 
###

invisible(tryCatch({get_libs()}, error = function(err) {
  end_pantera(paste("Installing libs error ->",geterrmessage()))
  }))
invisible(tryCatch({opt <- read_pars()}, error = function(err) {
  end_pantera(paste("Reading parameters error ->",geterrmessage()))
  }))
if (opt$log_file == FALSE) {
  lx <- function(x){return(NULL)}
}

# Internal globals 
zones_interval <- 60
zones_interval_overlap <- 30
current <- getwd()
segments_unique <- data.table(name = as.character(),
                              seq = as.character())

invisible(tryCatch({init_pantera()}, error = function(err) {
  end_pantera(paste("Initialization error ->",geterrmessage()))
  }))
invisible(tryCatch({segments_unique <- read_poly()}, error = function(err) {
  end_pantera(paste("GFA files missing ->",geterrmessage()))
  }))

setwd(opt$output_folder)

invisible(tryCatch({find_repeats()}, error = function(err) {
  end_pantera(paste("Error round 1 ->",geterrmessage()))
  }))
invisible(tryCatch({cluster_results()}, error = function(err) {
  end_pantera(paste("Error round 2 ->",geterrmessage()))
  }))
invisible(tryCatch({final <- classify_tes()}, error = function(err) {
  end_pantera(paste("Error Classification ->",geterrmessage()))
  }))
invisible(tryCatch({stats_tes()}, error = function(err) {
  end_pantera(paste("Error Stats ->",geterrmessage()))
  }))
if (opt$debug == FALSE) {
invisible(tryCatch({cleanup()}, error = function(err) {
  end_pantera(paste("Error Cleanning up ->",geterrmessage()))
  }))
}


