library(getopt)
library(methods)
suppressMessages(library(stringr))
suppressMessages(library(data.table))
suppressMessages(library(Biostrings))
suppressMessages(library(sybil))
library(htmltab)
options(error=traceback)

# get options
spec <- matrix(c(
  'input.dir', 'i', 1, "character", "Folder containing predicted traits.",
  'help' , 'h', 0, "logical", "help"
), ncol = 5, byrow = T)

opt <- getopt(spec)

# Help Screen
if ( !is.null(opt$help) | is.null(opt$input.dir) ){
  cat(getopt(spec, usage=TRUE))
  
  cat("\n")
  cat("Details:\n")
  q(status=1)
}

input.dir <- opt$input.dir
if( !dir.exists(input.dir) ) stop("Input directory not found.")

avail.dir <- list.dirs(input.dir, recursive=F, full.names=F)
cat("found:", avail.dir, "in:", input.dir, "\n")


# Virulence factors, resistance genes (abricate)
if("abricate" %in% avail.dir){
    files.abricate.vfdb <- list.files(paste0(input.dir, "/abricate"), pattern="*_vfdb.tbl", full.names=T)
    abricate.vfdb.dt <- rbindlist( Map(cbind, lapply(files.abricate.vfdb, data.table::fread), org = str_remove_all(basename(files.abricate.vfdb), "abricate_|_vfdb.tbl")) )
    files.abricate.resfinder <- list.files(paste0(input.dir, "/abricate"), pattern="*_resfinder.tbl", full.names=T)
    abricate.resfinder.dt <- rbindlist( Map(cbind, lapply(files.abricate.resfinder, data.table::fread), org = str_remove_all(basename(files.abricate.resfinder), "abricate_|_resfinder.tbl")) )
    write.csv(abricate.vfdb.dt, gzfile(paste0(input.dir, "/abricate_vfdb.csv.gz")), row.names = FALSE)
    write.csv(abricate.resfinder.dt, gzfile(paste0(input.dir, "/abricate_resfinder.csv.gz")), row.names = FALSE)

    # write.csv(abricate.resfinder.dt, paste0(input.dir, "/abricate_resfinder.csv.gz"))
}

# Natural products (antismash)
if("antismash" %in% avail.dir){
    files.antismash <- list.files(paste0(input.dir, "/antismash"), pattern="index.html", full.names=T,recursive=T)	
    antismash.dt <- data.table()
    for(f in files.antismash){
        no.result <- any(grepl("No results found on input", readLines(f, warn=F)))
        if(no.result) next
        new.dt <- tryCatch(
        expr = { return(data.table(htmltab(f,which="//div[@id='compact-record-table']",rm_nodata_cols=F), org = str_replace(basename(f), "_index.html$", ""))) },
        error = function(e){ return(data.table(htmltab(f,which="//div[@class='overview-layout']", rm_nodata_cols=F),org = str_replace(basename(f), "_index.html$", ""))) },
            finally = { })
        colnames(new.dt) <- c("Region","Type","From","To","Most.similar.known.cluster","Most.similar.known.cluster.1","Similarity","org")
        antismash.dt <- rbind(antismash.dt, new.dt)
    }
    write.csv(antismash.dt, gzfile(paste0(input.dir, "/antismash.csv.gz")), row.names = FALSE)
    # write.csv(antismash.dt, paste0(input.dir, "/antismash.csv.gz"))
}


# 16S genes (barrnap)
if("barrnap" %in% avail.dir){
    files.barrnap <- list.files(paste0(input.dir, "/barrnap"), pattern="*_16S.fna", full.names=T)
    barrnap.lst <- lapply(files.barrnap, readDNAStringSet)
    barrnap.dt <- data.table(org=str_remove(basename(files.barrnap),"_16S.fna"), rrna.copy=sapply(barrnap.lst, function(x){length(x[grep("16S",names(x))])}))
    fwrite(barrnap.dt, paste0(input.dir, "/barrnap.csv.gz"))
}

# Carbohydrate active enzymes (dbcan)
if ("dbcan" %in% avail.dir) {
    files.dbcan <- list.files(paste0(input.dir, "/dbcan"), pattern="overview.txt", full.names=TRUE, recursive=TRUE)
    dbcan.dt <- rbindlist(
        Map(
            function(dt, file) {
                dt[, org := tools::file_path_sans_ext(basename(file))] 
            dt
            },
            lapply(files.dbcan, data.table::fread),
            files.dbcan
        )
    )
    files.dbcan.sub <- list.files(paste0(input.dir, "/dbcan"), pattern="dbcan-sub.hmm.out", full.names=TRUE, recursive=TRUE)
    dbcan.sub.dt <- rbindlist(
        Map(
            function(dt, file) {
                dt[, org := tools::file_path_sans_ext(basename(file))] 
            dt
            },
            lapply(files.dbcan, data.table::fread),
            files.dbcan
        )
    )
    write.csv(dbcan.dt, gzfile(paste0(input.dir, "/dbcan.csv.gz")), row.names=FALSE)
    write.csv(dbcan.sub.dt, gzfile(paste0(input.dir, "/dbcan_sub.csv.gz")), row.names=FALSE)
}


# cluster of orthologous groups (eggnog-mapper)
if("eggnog" %in% avail.dir){
    files.eggnog <- list.files(paste0(input.dir, "/eggnog"), pattern="emapper.annotations", full.names=T,recursive=T)
    eggnog.dt <- rbindlist( Map(cbind, lapply(files.eggnog, function(f) data.table::fread(cmd=paste("grep -v '^##'", f))), org = str_extract(files.eggnog,"(?<=eggnog/).*(?=.emapper)")) )
    fwrite(eggnog.dt, paste0(input.dir, "/eggnog.csv.gz"))
}

# Gut gene cluster (gutsmash)
if("gutsmash" %in% avail.dir){
    files.gutsmash <- list.files(paste0(input.dir, "/gutsmash"), pattern="index.html", full.names=T,recursive=T)
    no.result.idx <- sapply(files.gutsmash, function(f) any(grepl("No results found on input", readLines(f, warn=F))))
    gutsmash.dt <- rbindlist( Map(cbind, lapply(files.gutsmash[!no.result.idx], htmltab,which=1,rm_nodata_cols=F),org = str_replace(basename(files.gutsmash[!no.result.idx]), "_index.html$", "")) )
    write.csv(gutsmash.dt, gzfile(paste0(input.dir, "/gutsmash.csv.gz")), row.names = FALSE)
}

                            
# Kegg (kofam)
if("kofam" %in% avail.dir){
    files.kofam <- list.files(paste0(input.dir, "/kofam"), pattern="*.txt", full.names=T,recursive=T)
    kofam.dt <- rbindlist( Map(cbind, data.table(t(lapply(files.kofam, readLines))), org = str_remove(basename(files.kofam),".txt")) )
    kofam.dt[,c("gene","kegg") :=tstrsplit(V1,"\t", fixed=F)]
    setnames(kofam.dt, old="V2", new="org"); kofam.dt[,V1:=NULL]
    fwrite(kofam.dt, paste0(input.dir, "/kofam.csv.gz"))
}

# metabolic network & model (gapseq)
if("gapseq" %in% avail.dir){
    files.gapseq.pwy <- list.files(paste0(input.dir, "/gapseq"), pattern="*-Pathways.tbl", full.names=T,recursive=T)
    gapseq.pwy.dt <- rbindlist( Map(cbind, lapply(files.gapseq.pwy, data.table::fread), org = str_extract(basename(files.gapseq.pwy),".*(?=-all-Pathways.tbl?)")) )
    files.gapseq.med <- list.files(paste0(input.dir, "/gapseq"), pattern="*-medium.csv", full.names=T,recursive=T)
    gapseq.med.dt <- rbindlist( Map(cbind, lapply(files.gapseq.med, data.table::fread), org = str_extract(basename(files.gapseq.med),".*(?=-medium.csv)")) )
    files.gapseq.cs <- list.files(paste0(input.dir, "/gapseq"), pattern="*-cs.tbl", full.names=T,recursive=T)
    gapseq.cs.dt <- rbindlist( Map(cbind, lapply(files.gapseq.cs, data.table::fread), org = str_extract(basename(files.gapseq.cs),".*(?=-cs.tbl)")) )
    files.gapseq.ferm <- list.files(paste0(input.dir, "/gapseq"), pattern="*-ferm.tbl", full.names=T,recursive=T)
    gapseq.ferm.dt <- rbindlist( Map(cbind, lapply(files.gapseq.ferm, data.table::fread), org = str_extract(basename(files.gapseq.ferm),".*(?=-ferm.tbl)")) )
    # Models
    files.gapseq.mod <- list.files(paste0(input.dir, "/gapseq"), pattern = "*.RDS", full.names = TRUE, recursive = TRUE)
    gapseq.mod.lst <- lapply(files.gapseq.mod, readRDS)
    names(gapseq.mod.lst) <- sapply(gapseq.mod.lst, function(x) {
        if ("mod_id" %in% slotNames(x)) x@mod_id else basename(files.gapseq.mod)[[1]]
    })
    write.csv(gapseq.pwy.dt, gzfile(paste0(input.dir, "/gapseq_pwy.csv.gz")), row.names = FALSE)
    write.csv(gapseq.med.dt, gzfile(paste0(input.dir, "/gapseq_med.csv.gz")), row.names = FALSE)
    write.csv(gapseq.cs.dt, gzfile(paste0(input.dir, "/gapseq_cs.csv.gz")), row.names = FALSE)
    write.csv(gapseq.ferm.dt, gzfile(paste0(input.dir, "/gapseq_ferm.csv.gz")), row.names = FALSE)
    saveRDS(gapseq.mod.lst, paste0(input.dir, "/gapseq_models.RDS"), compress = "xz")
}


if("bakta" %in% avail.dir){
    files.bakta <- list.files(paste0(input.dir, "/bakta"), pattern="*.txt", full.names=T,recursive=T)
    bakta.dt <- rbindlist( Map(cbind, lapply(files.bakta, data.table::fread, sep=":", fill=T), org = str_remove(basename(files.bakta),".txt")) )
    files.bakta.ext <- grep("hypothetical", list.files(paste0(input.dir, "/bakta"), pattern="*.tsv", full.names=T,recursive=T), invert=T, value=T)
    bakta.ext.dt <- rbindlist( Map(cbind, lapply(files.bakta.ext, data.table::fread), org = str_remove(basename(files.bakta.ext),".tsv")) )
    # amino acid sequences are needed for grodon (codon usage bias)
    files.bakta.prot <- grep("hypothetical", list.files(paste0(input.dir, "/bakta"), pattern="*.faa", full.names=T,recursive=T), value=T, invert=T)
    bakta.prot.seq <- lapply(files.bakta.prot, readAAStringSet)
    bakta.prot.dt <- rbindlist( Map(cbind, data.table(t(lapply(bakta.prot.seq, names))), org = str_remove(basename(files.bakta.prot),".faa")) )
    setnames(bakta.prot.dt, old=c("V1","V2"), new=c("protein","org"))
    fwrite(bakta.dt, paste0(input.dir, "/bakta.csv.gz"))
    fwrite(bakta.ext.dt, paste0(input.dir, "/bakta_ext.csv.gz"))
    fwrite(bakta.prot.dt, paste0(input.dir, "/bakta_prot.csv.gz"))
}

if("grodon" %in% avail.dir){
    files.grodon <- list.files(paste0(input.dir, "/grodon"), pattern="*.csv", full.names=T,recursive=F)
    grodon.dt <- rbindlist( Map(cbind, lapply(files.grodon, data.table::fread), org = str_remove(basename(files.grodon),".csv")) )
    fwrite(grodon.dt, paste0(input.dir, "/grodon.csv.gz"))
}

if("platon" %in% avail.dir){
    files.platon <- list.files(paste0(input.dir, "/platon"), pattern="*.tsv", full.names=T,recursive=F)
    if(length(files.platon)>0){
        platon.dt <- rbindlist( Map(cbind, lapply(files.platon, data.table::fread), org = str_remove(basename(files.platon),".tsv")) )[!is.na(ID)]
    }else platon.dt <- data.table()
    fwrite(platon.dt, paste0(input.dir, "/platon.csv.gz"))
}
