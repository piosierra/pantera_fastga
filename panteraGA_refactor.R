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
# 0.7.0 Numerous improvements in filters and Unknown recovery

# TODO 
# Create an special type of TSD check for satellites, when the TSD maps to the
# other edge (as series of tandem repeats)
# 
# Remove also "identical" sequences, not only fragments. 
#

pantera_version <- "0.6.2"
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
    install.packages("this.path") #
  }
  if (suppressPackageStartupMessages(!require("getopt", quietly = TRUE))) {
    lx("Intalling package [getopt]")
    install.packages("getopt") #
  }
  if (suppressPackageStartupMessages(!require("parallel", quietly = TRUE))) {
    lx("Intalling package [parallel]")
    install.packages("parallel") #
  }
  if (suppressPackageStartupMessages(!require("ips", quietly = TRUE))) {
    lx("Intalling package [ips]")
    install.packages("ips") #
  }
  if (suppressPackageStartupMessages(!require("BiocManager", quietly = TRUE))) {
    lx("Intalling package [BiocManager]")
    install.packages("BiocManager") #
  }
  if (suppressPackageStartupMessages(!require("Biostrings", quietly = TRUE))) {
    lx("Intalling package [Biostrings]")
    BiocManager::install("Biostrings") #
  }
  if (suppressPackageStartupMessages(!require("DECIPHER", quietly = TRUE))) {
    lx("Intalling package [DECIPHER]")
    BiocManager::install("DECIPHER") #
  } 
  if (suppressPackageStartupMessages(!require("stringi", quietly = TRUE))) {
    lx("Intalling package [stringi]")
    install.packages("stringi") #
  }
  if (suppressPackageStartupMessages(!require("ape", quietly = TRUE))) {
    lx("Intalling package [ape]")
    install.packages("ape") #
  }
  if (suppressPackageStartupMessages(!require("LncFinder", quietly = TRUE))) {
    lx("Intalling package [LncFinder]")
    install.packages("LncFinder") #
  }
  if (suppressPackageStartupMessages(!require("xgboost", quietly = TRUE))) {
    lx("Intalling package [xgboost]")
    install.packages("xgboost") #
  }
  if (suppressPackageStartupMessages(!require("quantmod", quietly = TRUE))) {
    lx("Intalling package [quantmod]")
    install.packages("quantmod") #
  }
  if (suppressPackageStartupMessages(!require("data.table", quietly = TRUE))) {
    lx("Intalling package [data.table]")
    install.packages("data.table") #
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
      "flanking",       "f", 1, "integer",   # Min length of flanking sequences [automatic]
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
    opt$min_size <- 200
  }
  
  if (is.null(opt$max_size)) {
    opt$max_size <- 20000
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
    # Length of polyA.
    opt$pas <- 10
  }
  
  if (is.null(opt$pag)) {
    # Max gap in polyA.
    opt$pag <- 5
  }
  
  if (is.null(opt$cl_size)) {
    # Maximum size of initial cluster to process
    opt$cl_size <- 300
  }
  
  if (is.null(opt$flanking)) {
    # Min length of mapping flanking sequences to define an SV. 100 performs automatic adjustment
    opt$flanking <- 100
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

# Writes a fasta file
wfasta <- function(fasta_data, f) {
  fileconn <- file(f)
  writeLines(t(fasta_data), fileconn)
  close(fileconn)
}

## Code taken from https://github.com/etam4260/kneedle/blob/main/R/kneedle.R
## Give atribution!
kneedle <- function(x, y, decreasing, concave, sensitivity = 1) {
  # Make sure inputs are correct
  if(length(x) == 0 || length(y) == 0) {
    stop("Make sure size of both inputs x and y are greater than 0")
  }
  if(typeof(x) == "list"|| typeof(y) == "list" || is.data.frame(x) ||
     is.data.frame(y) || is.array(x) || is.array(y) || is.matrix(x) || is.matrix(y)) {
    stop("Make sure both inputs x and y are vectors")
  }
  if(length(x) != length(y)) {
    stop("Make sure size of both inputs x and y are equal")
  }
  
  data <- matrix(unlist(list(x, y)), ncol = 2)
  # This decreasing value has nothing to do with the user inputted value.
  data <- data[order(data[,1], decreasing = FALSE), ]
  
  
  # If both decreasing and concave are not specified, then this algorithm will
  # take a guess at those parameters instead of defaulting to certain values.
  # One method is to take the derivative of the starting point to the ending point
  # from min x value to max x value. This should be similar to taking the average of all
  # derivatives from xi to xi+1 where i = 1 to i = length(xvalues)
  if(missing(decreasing)) {
    # Increasing discrete data
    if( (data[(nrow(data)), 2] - data[1, 2]) >= 0 ) {
      decreasing = FALSE
      # Decreasing discrete data
    } else {
      decreasing = TRUE
    }
  }
  
  # To determine concavity we need to look at the second derivative of the
  # entire set of discrete data from xi to xi+1 where i = 1 to i = length(xvalues)
  # Taking the average of all the second derivatives, if greater or equal to 0
  # then concave up. If less than 0 then concave down.
  if(missing(concave)) {
    secondderiv <- diff(diff(data[, 2]) / diff(data[, 1]))
    if(mean(secondderiv) > 0) {
      concave = TRUE
    } else {
      concave = FALSE
    }
  }
  
  maxy <- max(y)
  miny <- min(y)
  maxx <- max(x)
  minx <- min(x)
  data[ ,1] <- (data[, 1]- min(data[, 1]))/(max(data[ ,1])- min(data[, 1]))
  data[ ,2] <- (data[, 2]- min(data[, 2]))/(max(data[ ,2])- min(data[, 2]))
  
  if(concave && !decreasing) {
    differ <- abs(c(data[ ,2] - data[ ,1]))
  } else if(concave && decreasing) {
    differ <- abs(c(data[ ,2] - (1 - data[ ,1])))
  } else if(!concave && !decreasing) {
    differ <- abs(c(data[ ,2] - data[ ,1]))
  } else if(!concave && decreasing) {
    differ <- abs(c(data[ ,2] - (1 - data[ ,1])))
  }
  
  
  peak.indices <- findPeaks(differ) - 1
  
  data <- cbind(data, differ)
  
  diffx = diff(data[, 1])
  T.lm.x.s <- sensitivity * mean(diffx)
  knee = NULL
  
  for(i in 1:length(peak.indices)) {
    T <- data[peak.indices[i] ,3] - (T.lm.x.s)
    
    y.value <- data[peak.indices[i] ,3]
    
    for(j in (peak.indices[i]):if(i+1 < length(peak.indices)) peak.indices[i+1] else length(differ)) {
      if(differ[j] < T) {
        knee = peak.indices[i];
        break;
      }
    }
    if(!is.null(knee)) {
      break;
    }
  }
  
  # Returns the x,y coordinate values
  x <- ((maxx - minx) * (data[knee, 1])) + minx
  y <- ((maxy - miny) * (data[knee, 2])) + miny
  
  return(c(as.numeric(x),as.numeric(y)))
}

# Reads the two ONEcode files extracted from a 1aln file with svfind
parseONEview <- function(f) {
  lx(paste("Procesing file =", f))
  check <- file.size(f)
  if (check == 0 | is.na(check)) {
    end_pantera(".1aln file missing")
  } else {
   system(paste0("svfind -x 12 -f ", opt$flanking, " -a ",opt$output_folder,"/",opt$lib_name,"hapa -b ",opt$output_folder,"/",opt$lib_name,"hapb ", f, " > /dev/null 2>&1"))
    ova <- fread(cmd = paste0("ONEview ",opt$output_folder,"/",opt$lib_name,"hapa"), fill = T)
    ovb <- fread(cmd = paste0("ONEview ",opt$output_folder,"/",opt$lib_name,"hapb"), fill = T)
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
    if (opt$flanking == 100) {
    kneedle_data <- data.table(flanks=c(ovall$flank_r,ovall$flank_l))[,b:=flanks %/% 100][,.N,b][order(b)]
    opt$flanking <- tryCatch({kneedle(kneedle_data$b,kneedle_data$N)[1]*100},error = function(msg){return(2001)})
    }
    lx(paste("Flanking sequences cutpoint = ", opt$flanking))
    ovcut <- ovall[flank_l>=opt$flanking & flank_r>=opt$flanking]
    ### Remove all but one copy when the same segment can map to several places.
    ovall <-ovall[order(file,v_name,v_start)]
    ovall[,ov:=shift(v_start,-1), by=v_name]
    ovall[,ov:=ov-v_end]
    ovall[is.na(ov) ,ov:=1]
    ovall[ov < -seq_len ,ov:=1]
    ovall <- ovall[ov>0]
   
    unlink(paste0(opt$output_folder,"/",opt$lib_name,"hapa"))
    unlink(paste0(opt$output_folder,"/",opt$lib_name,"hapb"))
    
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
lx(paste("1aln file:", opt$sv_file))
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
          lx(paste0("Num. of segments zone ",start," - ",end," : ", nrow(segment_u)))
          segment_sets <- split(segment_u, 
                                rep(1:(nrow(segment_u) %/% opt$cl_size +1), 
                                length.out = nrow(segment_u), 
                                each = ceiling(nrow(segment_u) / 
                                               (nrow(segment_u) %/% 
                                                opt$cl_size +1))))
          candidates <- data.table(name = as.character(), 
                                   seq = as.character())
          for (ss in 1:length(segment_sets)) {
       #     lx(paste0("Num. of segments zone ", start," - ", end," : cluster : ",ss))
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
          
       #     lx(paste("# of clusters ", start, "-", end, ":", length(unique(good_clust) ))) 
         
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
 #     lx(paste("Procesing segments:", start, "-", end))
      if (nrow(segments_unique[len >= start][len <= end])> 0) {
        segment_u <- segments_unique[len >= start][len <= end]
        segment_u <- segment_u[!duplicated(segment_u)]
  #      lx(paste("Number of sequences:", nrow(segment_u)))
          if (nrow(segment_u) >= opt$min_cl)  {
            segment_sets <- split(segment_u, 
                                  rep(1:(nrow(segment_u) %/% opt$cl_size +1), 
                                      length.out = nrow(segment_u), 
                                      each = ceiling(nrow(segment_u) / 
                                                     (nrow(segment_u) %/%
                                                      opt$cl_size +1))))
    #        lx(paste("Segment sets", start, "-", end, ":", length(segment_sets)))
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
              seqs <- ape::as.DNAbin(clust_temp)
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
                cons <- toupper(Biostrings::consensusString(Biostrings::DNAStringSet(alim), 
                                                ambiguityMap = "-", 
                                                threshold = cons_threshold))
                
                ### Exploring how to find the TSDs. The cons gives as already the positions. Maybe get pos +-1,checkt that one, and then extend to both sides?
                
                cMatrix <- Biostrings::consensusMatrix(Biostrings::DNAStringSet(alim))
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
              #  tsds <- unlist(lapply(2:13,function(x){sum(diag(adist(substr(edges$r,1,x),stri_reverse(substr(edges$l,1,x)))))}))
                tsds2 <- unlist(lapply(2:13,function(x){sum(substr(edges$r,1,x) == stri_reverse(substr(edges$l,1,x)))}))
                tsd_len <- min(which(tsds2==max(tsds2),tsds2))+1
                tsds_conf <- round(tsds2[tsd_len-1]/nrow(edges),2)
                tsds_motif <- Biostrings::consensusString(Biostrings::DNAStringSet(c(substr(edges$r,1,tsd_len),stri_reverse(substr(edges$l,1,tsd_len)))), ambiguityMap=IUPAC_CODE_MAP,
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
  opt$cl_size <- 1000
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
      " --outseq temp.orfs -minsize 300 &>/dev/null; blastp -num_threads ",opt$threads," -query temp.orfs -db ",
      this.dir(),
      "/libs/RepeatPeps.lib -outfmt 6 -evalue 1e-1 > orfs.tbl"
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
  final[Probability < 0.6, name := paste0(gsub("#.*","",name), "#Unknown", collapse = ""), by=.I]
  final <- final[!duplicated(final$seq)]
  fwrite(final, paste0(file, ".statspre"), sep = "\t")
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
                               " -evalue 500 -outfmt 6 -word_size 11 -gapopen 4 -gapextend 1 -reward 1 -penalty -1"), header = F)
  if (nrow(te_data) > 0) {
    colnames(te_data) <-c("qseqid","sseqid", "pident" , "length","mismatch", 
                          "gapopen","qstart","qend","sstart","send", 
                          "evalue", "bitscore")
    te_data_mix <- te_data[qseqid != sseqid] # Non self matches, to deal with later
    te_data_rep <- te_data[qseqid == sseqid & length >20 & 
                             (qstart != sstart | qend != send)][,.(N=.N,maxRep=max(length)),qseqid] # large self matches that are not trivial. To filter with MaxRep.
    te_data  <- te_data[!(qseqid == sseqid & qstart == sstart & qend == send)] # Self matches that are not the whole element, to check for TIR and LTR
    ### TIR, LTR detection
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
    te_data$type <- "LTR"
    te_data[check_type==T,type:= "TIR"]
    ### TIR, LTR detection END
    
    # Mark as PASS LTR and TIR elements matching their class.
    good_ltr <- te_data[grepl("LTR",qseqid) & type == "LTR" & lgap < 8 & rgap < 8 & length > 100]$qseqid
    good_tir <- te_data[grepl("DNA",qseqid) & type == "TIR" & lgap < 8 & rgap < 8]$qseqid
    

  

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
   
   good_line <- tes[grepl("LINE",name) & orf1 > 1600 & (pa < 10 | is.na(type) | (lgap >10 & rgap > 10)) ]$name
   
   # Find elements which share TIR or LTR to good ones.
   te_data_mix <- merge(te_data_mix,
                        tes[,c("name","lente", "orf1", "orf2", "orf3")]
                        [,name:=gsub(">","",name)], by.x = "qseqid",
                        by.y="name")
   te_data_mix <- merge(te_data_mix,
                        tes[,c("name","lente", "orf1", "orf2", "orf3")]
                        [,name:=gsub(">","",name)], by.x = "sseqid",
                        by.y="name")
   te_data_mix[,lqgap:=.(min(qstart,qend)-1), by=.I]
   te_data_mix[,rqgap:=.(lente.x-max(qstart,qend)), by=.I]
   te_data_mix[,lsgap:=.(min(sstart,send)-1), by=.I]
   te_data_mix[,rsgap:=.(lente.y-max(sstart,send)), by=.I]
   ## Remove fully contained in good ones
   te_data_mix[,cover:=length/min(lente.x,lente.y), by=.I]
   te_data_mix[,small:= (min(lente.x,lente.y) / max(lente.x,lente.y)) < 0.5,by=.I]
   smaller <- te_data_mix[cover>0.98 & pident>98 & small == T]
   smaller_list <- unique(c(smaller[orf1.y>=orf1.x]$qseqid,smaller[orf1.x>orf1.y]$sseqid))
   tes <- tes[!(name %in% smaller_list)]
   te_data_mix <- te_data_mix[!(sseqid %in% smaller_list) & !(qseqid %in% smaller_list)]
   lx(paste("Removed fragments:", length(smaller_list)))
   
  # te_data_mix <- te_data_mix[!grepl("Unknown", sseqid)]
   te_data_mix <- te_data_mix[lsgap<10 | rsgap<10]
   te_data_mix <- te_data_mix[lqgap<10 | rqgap<10]
   
   ## Generate list of LINEs that are subsequence of a better one
   te_data_mix_line <- te_data_mix[grepl("LINE",qseqid) & grepl("LINE",sseqid)]
   te_data_mix_line[,large:=qseqid]
   te_data_mix_line[lente.y>lente.x,large:=sseqid]
  
   te_data_mix_line <- te_data_mix_line[cover >0.9]
   small_line <- c()
   for (l in tes[name %in% good_line][order(-lente)]$name) {
     matches <- unique(c(te_data_mix_line[large==l]$qseqid,te_data_mix_line[large==l]$sseqid))
     matches <- matches[matches!=l]
     small_line <- c(small_line,matches)
     te_data_mix_line <- te_data_mix_line[!(qseqid %in% matches | sseqid %in% matches)]
   }
   
   ### Find TIR elements that can be reclassified
   te_data_mix_tir <- te_data_mix[(qseqid %in% good_tir | sseqid %in% good_tir) & lsgap < 8 & rsgap < 8 & lqgap < 8 & rqgap <8 ]
   # tir_reco1 <- te_data_mix_tir[grepl("Unknown",qseqid)]
   # tir_reco2 <- te_data_mix_tir[grepl("Unknown",sseqid)]
    tir_reco1 <- te_data_mix_tir[!grepl("DNA",qseqid) & grepl("DNA",sseqid)]
    tir_reco2 <- te_data_mix_tir[!grepl("DNA",sseqid) & grepl("DNA",qseqid)]
   tir_reco <- data.table(name=c(tir_reco1$qseqid,tir_reco2$sseqid), sf=gsub(".*#","",c(tir_reco1$sseqid,tir_reco2$qseqid)))
   tir_reco <- tir_reco[!duplicated(tir_reco)]
   lx(paste("Unknown elements reclassified as DNA:", nrow(tir_reco)))
   tir_reco[,new := paste0(gsub("#.*","",name),"#",sf,collapse=""), by=.I]
   tir_merge <- merge(tes[,1],tir_reco,all=T)
   tir_merge[!is.na(new),name:=new]
   tes$name <- tir_merge$name
   
   ### Find LTR elements that can be reclassified
   te_data_mix_ltr <- te_data_mix[(qseqid %in% good_ltr | sseqid %in% good_ltr) & ((lsgap < 8  & rsgap < 8) | (lqgap < 8 & rqgap <8)) ]
   ltr_reco1 <- te_data_mix_ltr[grepl("Unknown",qseqid)]
   ltr_reco2 <- te_data_mix_ltr[grepl("Unknown",sseqid)]
   ltr_reco <- data.table(name=c(ltr_reco1$qseqid,ltr_reco2$sseqid), sf=gsub(".*#","",c(ltr_reco1$sseqid,ltr_reco2$qseqid)))
   ltr_reco <- ltr_reco[!duplicated(ltr_reco)]
   lx(paste("Unknown elements reclassified as LTR:", nrow(ltr_reco)))
   ltr_reco[,new := paste0(gsub("#.*","",name),"#",sf,collapse=""), by=.I]
   ltr_merge <- merge(tes[,1],ltr_reco,all=T)
   ltr_merge[!is.na(new),name:=new]
   tes$name <- ltr_merge$name
   
   tes[,pass:=F]
   tes[,sf:=gsub(".*#","",name)]
   tes[name %in% good_line[!(good_line %in% small_line)], pass:=T]
   tes[grepl("#DIRS",name) & orf1 > 2000 & lente < 10000, pass:=T]
   tes[grepl("Crypton",name), pass:=T]
   tes[grepl("#PLE",name), pass:=T]
   tes[grepl("#SINE",name) & lente < 450, pass:=T]
   tes[name %in% good_ltr, pass:=T]
   tes[name %in% ltr_reco$new, pass:=T]
   tes[name %in% good_tir, pass:=T]
   tes[name %in% tir_reco$new, pass:=T]
   tes[grepl("#RC",name) & orf1> 2000, pass:=T]
   tes[cluster>5 , pass:=T]
   
   # This should go at the end of the filters
   tes[is.na(pa),pa:=NA]
  # tes[,TR1:=maxRep>length & maxRep>1000]
   tes[,TR2:=maxRep>lente/2.5]
  # tes[is.na(TR1), TR1:=F]
   tes[is.na(TR2), TR2:=F]
  # tes[TR1!=T & TR2!=T, pass:=T]
   tes[TR2 == T & type == "LTR" & maxRep > 100, pass == F]
   tes[type == "TIR" & lgap < 3 & rgap< 3, pass:=T]
   
   # Reclassification by TSDs
   tes[grepl("#Unknown",name) & type == "TIR" & tsd_l == 8 & tsd_c >0.3 & lgap < 8 & rgap < 8,  `:=`(name=paste0(gsub("#.*","",name),"#DNA/hAT",collapse=""),pass=T), by=.I]
   tes[grepl("#Unknown",name) & type == "TIR" & lgap < 8 & rgap < 8,  `:=`(name=paste0(gsub("#.*","",name),"#DNA",collapse=""),pass=T), by=.I]
   tes[grepl("#Unknown",name) & tsd_l == 4 & tsd_c >0.3 & type != "TIR", `:=`(name=paste0(gsub("#.*","",name),"#LTR",collapse=""),pass=T), by=.I]
   tes[grepl("#Unknown",name) & tsd_l == 5 & tsd_c >0.3 & type != "TIR", `:=`(name=paste0(gsub("#.*","",name),"#LTR",collapse=""),pass=T), by=.I]
   tes[grepl("#Unknown",name) & tsd_l == 6 & tsd_c >0.3 & type != "TIR", `:=`(name=paste0(gsub("#.*","",name),"#LTR",collapse=""),pass=T), by=.I]
   tes[grepl("#Unknown",name) & tsd_l > 3 & tsd_l < 7 & tsd_c >0.3 & type == "LTR" & length > 150 & lgap < 8 & rgap < 8 & TR2 == F, `:=`(name=paste0(gsub("#.*","",name),"#LTR",collapse=""),pass=T), by=.I]
   
   # SINE reclassification by pA and size
   tes[grepl("#Unknown",name) & lente < 450 & pa < 5, `:=`(name = paste0(gsub("#.*","",name),"#SINE",collapse=""),pass=T), by=.I]
   
 
  lx(paste("TRs discards:", nrow(tes[pass==F])))
  }
  setwd("..")
  unlink(temp, recursive = TRUE)
  tes <- tes[order(-lente)]
  tes[,name:=paste0(">",name, collapse = ""), .I]
  wfasta(tes[, c("name", "seq")], paste0(opt$lib_name, "-pantera-final.fa"))
  wfasta(tes[pass == T, c("name", "seq")], paste0(opt$lib_name, "-pantera-pass.fa"))
  wfasta(tes[pass == F, c("name", "seq")], paste0(opt$lib_name, "-pantera-discards.fa"))
  stats_data <- tes[,c("name", "cluster", "tsd_l","tsd_c","tsd_m","pass","lente","type","pa",
                       "length","lgap","rgap","orf1","orf2","orf3","maxRep")]
  ### Last fixes
  stats_data[tsd_m == "" | is.na(tsd_m), tsd_c :=0]
  
  colnames(stats_data) <- c("name", "cluster_size", "TSD_length","TSD_confidence","TSD_motif","pass","TE_len", 
                            "struct_type","polyA","struct_len","left_gap",
                            "right_gap","orf1","orf2","orf3","maxTR")

  fwrite(stats_data, paste0(opt$lib_name, "-pantera-final.stats"), 
         quote =  F, row.names = F, sep ="\t")
}

# Delete temporary files
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
  end_pantera(paste("Error reading 1aln ->",geterrmessage()))
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


