#!/usr/bin/env Rscript

# new version of pantercheck

options(warn = 0)


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


if (suppressPackageStartupMessages(!require("data.table", quietly = TRUE))) {
  install.packages("data.table")
}
if (suppressPackageStartupMessages(!require("LncFinder", quietly = TRUE))) {
  install.packages("LncFinder")
}
if (suppressPackageStartupMessages(!require("parallel", quietly = TRUE))) {
  install.packages("parallel")
}
if (suppressPackageStartupMessages(!require("stringi", quietly = TRUE))) {
  install.packages("stringi")
}

### Function to read a fasta file

ffasta <- function(f) {
  check <- file.size(f)
  if (check == 0 | is.na(check)) {
    return(data.table(name = as.character(), seq = as.character()))
  } else {
    fa_raw <- fread(f,
                    header = F,
                    fill = T,
                    sep = "\n")
    fa_raw[, h := grepl(">", V1)]
    fa_oneline <- fa_raw[, .(paste0(V1, collapse = "")), by = rleid(h)]
    return(data.table(name = fa_oneline[rep(c(TRUE, FALSE), length = .N)]$V1, seq =
                        fa_oneline[rep(c(FALSE, TRUE), length = .N)]$V1))
  }
}

### Function to write a fasta file

wfasta <- function(fasta_data, f) {
  fileconn <- file(f)
  writeLines(t(fasta_data), fileconn)
  close(fileconn)
}




###
### Main functions 
###

# Init function

  f <- commandArgs(trailingOnly=TRUE)[1]
  tes <- ffasta(f)
  temp <- paste0("temp",make.names(gsub(".*/","",f)))
  dir.create(temp)
  setwd(temp)
  tes[,lente:=nchar(seq)]
  tes[,name:=gsub(" .*","",name)]
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
  te_data <- fread(cmd= paste0("blastn -query tmp-tes -db tmp-tes -task blastn -num_threads 8 -evalue 500 -outfmt 6 -word_size 7 -gapopen 4 -gapextend 1 -reward 1 -penalty -1"), header = F)
  if (nrow(te_data) > 0) {
    colnames(te_data) <-c("qseqid","sseqid", "pident" , "length","mismatch", 
                          "gapopen","qstart","qend","sstart","send", 
                          "evalue", "bitscore")
    te_data_mix <- te_data[qseqid != sseqid]
    te_data_rep <- te_data[qseqid == sseqid & 
                             (qstart != sstart | qend != send)][,.N,qseqid]
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
    te_data[,pass:=T]
   #  te_data[grepl("LTR",qseqid) & !grepl("DIRS",qseqid) & 
  #          (orf1<1000 | lgap > 20 | rgap > 20 | type != "LTR" |
   #          lente < 3000 | lente > 16000), pass:=F]
    # te_data[grepl("LTR",qseqid) & grepl("DIRS",qseqid) & 
    #        (orf1<1000 | lgap > 20 | rgap > 20 | type != "TIR" | 
    #         lente < 3000 | lente > 16000), pass:=F]
   # te_data[grepl("DNA",qseqid) & (lgap > 20 | rgap > 20  | type != "TIR"), pass:=F]
   #  te_data[grepl("RC",qseqid) & (lente < 4000 |orf1<2000) , pass:=F]
   # te_data[grepl("LINE",qseqid) &  orf1 < 1600, pass:=F]
   # te_data[grepl("DIRS",qseqid) & (lente > 4000 | lente < 9000 | orf1 > 2000), pass:=T]
   # te_data[grepl("Crypton",qseqid), pass:=T]
    

    tes <- merge(tes[,name:=gsub(">","",name)], 
                 te_data[,c("qseqid","length","lgap","rgap", "type", "pass")], 
                          by.x = "name", by.y = "qseqid", all.x=T)
   # tes[is.na(pass), pass :=F]
    tes[is.na(type), rgap :=0]
    tes[is.na(type), lgap :=0]
    tes[is.na(type), length :=0]
    
    tes$pas <- unlist(lapply(stri_locate_all_regex(tes$seq,paste0(strrep("A",10),"|",strrep("T",10))),function(x) {min(unlist(x))-1}))
    tes$eas <- nchar(tes$seq)-unlist(lapply(stri_locate_all_regex(tes$seq,paste0(strrep("A",10),"|",strrep("T",10))),function(x) {max(unlist(x))}))
    tes[,pa:=min(pas,eas), by=.I]
    
    tes <- merge(tes,te_data_rep, by.x = "name", by.y = "qseqid", all.x=T)
    
    # LINE_check <-  te_data_mix[sseqid %in% te_data[pass==T & 
    #                            grepl("LINE",sseqid)]$sseqid &
    #                            !(qseqid %in% te_data[pass==T]$sseqid) &
    #                            ((length-lente.x)^2)<1000]
    # #  DNA_check <-  te_data_mix[sseqid %in% te_data[pass==T & 
    #                 # grepl("DNA",sseqid)]$sseqid & qseqid %in% te_data[pass==F 
    #                 # & type == "TIR"]$sseqid & (lqgap<20 | rqgap<20)]
    # #  DNA_check2 <-  te_data_mix[sseqid %in% te_data[pass==F & 
    #                 # type =="TIR"]$sseqid & qseqid %in% te_data[pass==T 
    #                 # & grepl("DNA",sseqid)]$sseqid ]
    # LTR_check <-  te_data_mix[sseqid %in% te_data[pass==T &
    #                           grepl("LTR",sseqid)]$sseqid &
    #                           !(qseqid %in% te_data[pass==T]$sseqid) &
    #                             ((length-lente.x)^2)<1000]
    # # TODO: Add search for NA LTR.
    # tes[(name %in% LINE_check$qseqid), pass:=F]
    # tes[(name %in% LTR_check$qseqid), pass:=F]
    # # tes[(name %in% te_data_rep[N>20]$qseqid & orf1 < 1000), pass:=F] # These should be parameters (was 100 and 600)
    # tes[(length > (lente / 2)), pass:=F]
   
  }
  setwd("..")
  unlink(temp, recursive = TRUE)
  tes <- tes[order(-pass,-lente)]
  tes_fa <- tes[, c("name", "seq")]
  stats_data <- tes[,c("name","lente", "pass","type", "pa",
                       "length","lgap","rgap","orf1","orf2","orf3","N")]
  colnames(stats_data) <- c("name","TE_len", "pass",
                            "struct_type", "polyA","struct_len","left_gap",
                            "right_gap","orf1","orf2","orf3","N")
  fwrite(stats_data, paste0(f, ".stats"),  quote =  F, row.names = F, sep ="\t")



