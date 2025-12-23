
library(getopt)
library(methods)
suppressMessages(library(stringr))
suppressMessages(library(stringi))
suppressMessages(library(data.table))
suppressMessages(library(cobrar))
suppressMessages(library(sybil))
options(error=traceback)
library(data.table)              #- added
library(R.utils)                  # -- added

# Define command-line options
spec <- matrix(c(
  'input.dir', 'i', 1, "character", "Folder containing predicted traits.",
  'db.path',   'd', 1, "character", "Path to gapseq dat folder.",
  'help',      'h', 0, "logical",   "Show help message and exit."
), ncol = 5, byrow = TRUE)

opt <- getopt(spec)

# Help Screen
if (!is.null(opt$help) | is.null(opt$input.dir) | is.null(opt$db.path)) {
  cat(getopt(spec, usage = TRUE))
  cat("\nDetails:\n")
  cat("This script reads fusion traits and loads pathway metadata from the specified gapseq database path.\n")
  q(status = 1)
}

# Validate input
input.dir <- opt$input.dir
db.path   <- opt$db.path

if (!dir.exists(input.dir)) stop("Input directory not found.")
if (!dir.exists(db.path))   stop("Gapseq dat/ folder not found.")

# Construct full paths
meta_file   <- file.path(db.path, "meta_pwy.tbl")
print(meta_file)
custom_file <- file.path(db.path, "custom_pwy.tbl")
print("opnedd!!")

# Validate files
if (!file.exists(meta_file))   stop("meta_pwy.tbl not found.")
if (!file.exists(custom_file)) stop("custom_pwy.tbl not found.")

# Load data
db.meta <- rbind(
  fread(meta_file),
  fread(custom_file)
)
# db.meta
abricate.vfdb.dt <- fread(paste0(input.dir, "/abricate_vfdb.csv.gz"))
abricate.resfinder.dt <- fread(paste0(input.dir, "/abricate_resfinder.csv.gz"))
antismash.dt <- fread(paste0(input.dir, "/antismash.csv.gz"))
barrnap.dt <- fread(paste0(input.dir, "/barrnap.csv.gz"))
dbcan.dt <- fread(paste0(input.dir, "/dbcan.csv.gz"))
dbcan.sub.dt <- fread(paste0(input.dir, "/dbcan_sub.csv.gz"))
eggnog.dt <- fread(paste0(input.dir, "/eggnog.csv.gz"))
gutsmash.dt <- fread(paste0(input.dir, "/gutsmash.csv.gz"))
kofam.dt <- fread(paste0(input.dir, "/kofam.csv.gz"))
gapseq.pwy.dt <- fread(paste0(input.dir, "/gapseq_pwy.csv.gz"))
gapseq.med.dt <- fread(paste0(input.dir, "/gapseq_med.csv.gz"))
gapseq.cs.dt <- fread(paste0(input.dir, "/gapseq_cs.csv.gz"))
gapseq.ferm.dt <- fread(paste0(input.dir, "/gapseq_ferm.csv.gz"))
gapseq.models <- readRDS(paste0(input.dir, "/gapseq_models.RDS"))
bakta.dt <- fread(paste0(input.dir, "/bakta.csv.gz"))
bakta.prot.dt <- fread(paste0(input.dir, "/bakta_prot.csv.gz"))
bakta.ext.dt <- fread(paste0(input.dir, "/bakta_ext.csv.gz"))
grodon.dt <- fread(paste0(input.dir, "/grodon.csv.gz"))
platon.dt <- fread(paste0(input.dir, "/platon.csv.gz"))

lht.dt <- data.table(org=grodon.dt$org)
print("db")
print(lht.dt)
# codon usage bias
lht.dt <- merge(lht.dt, grodon.dt[,.(GC,CUB,CPB,d), by="org"])
# Based upon d growth type is defined
lht.dt[,grodon_growth.type:=ifelse(d>5,"slow","fast")]
setnames(lht.dt,old=c("GC","CUB","CPB","d"),new=c("grodon_GC","grodon_CUB","grodon_CPB","grodon_d"))
print(lht.dt, "codon usage")
#write.csv(lht.dt, gzfile(paste0(input.dir, "GRODEN_MAGs.csv.gz")), row.names=FALSE)
# auxotrophies
aa.lst <-c("cpd00041","cpd00130","cpd00322","cpd00107","cpd00156","cpd00161","cpd00035","cpd00051","cpd00132","cpd00041","cpd00084","cpd00023","cpd00053","cpd00033","cpd00119","cpd00322","cpd00107","cpd00039","cpd00060","cpd00064","cpd00066","cpd00129","cpd00054","cpd00161","cpd00065","cpd00069","cpd00156","cpd00266")
vit.lst <- c("cpd01631","cpd00365","cpd03185","cpd00220","cpd00215","cpd03424","cpd00541","cpd00104","cpd00393","cpd00305","cpd00218","cpd00644","cpd04145","cpd01401")
lht.dt <- merge(lht.dt,  gapseq.med.dt[,list(gapseq_auxotrophy=sum(compounds%in%c(aa.lst,vit.lst)), gapseq_auxotrophy.aa=sum(compounds%in%aa.lst), gapseq_auxotrophy.vit=sum(compounds%in%vit.lst)), by=org], by="org")
print( "auxo")
print(lht.dt)
# annotation (bakta)
bakta.tmp <- dcast(data=bakta.dt[V1%in%c("Length","Count","GC","N50","N ratio","coding density","tRNAs","tmRNAs","rRNAs","ncRNAs","ncRNA regions", "CRISPR arrays","CDSs","pseudogenes","signal peptides","sORFs")], formula=org~V1, value.var="V2")
colnames(bakta.tmp) <- ifelse(colnames(bakta.tmp)=="org", colnames(bakta.tmp), paste0("bakta_", colnames(bakta.tmp)))
lht.dt <- merge(lht.dt, bakta.tmp, by="org")
print("bakta")
print(head(lht.dt))
# virulence
vfdb.categories.lst <- c("Adherence", "Antimicrobial activity/Competitive advantage", "Biofilm", "Effector delivery system", "Exotoxin", "Exoenzyme", "Immune modulation", "Invasion", "Motility", "Nutritional/Metabolic factor", "Regulation", "Stress survival", "Post-translational modification", "Others")
abricate.vfdb.dt[!is.na(GENE),category:=str_extract(PRODUCT, paste0(vfdb.categories.lst,collapse="|"))]
abricate.vfdb.tmp <- dcast(data=abricate.vfdb.dt, formula=org~category, fun=length)
for(col in setdiff(vfdb.categories.lst, colnames(abricate.vfdb.tmp))) abricate.vfdb.tmp[,(col):=0] # add missed categories with zero values
if("NA" %in% colnames(abricate.vfdb.tmp)) abricate.vfdb.tmp <- abricate.vfdb.tmp[,-"NA"]
colnames(abricate.vfdb.tmp) <- ifelse(colnames(abricate.vfdb.tmp)=="org","org", paste0("vfdb_",make.names(colnames(abricate.vfdb.tmp))))
lht.dt <- merge(lht.dt, abricate.vfdb.tmp, by="org")
print("virulence")
# print(lht.dt)
# resistance
abricate.resfinder.tmp <- abricate.resfinder.dt[,list(resfinder_hits=length(stri_remove_empty(GENE)), resfinder_genes=length(stri_remove_empty(unique(GENE))), resfinder_targets=length(stri_remove_empty(unique(unlist(str_split(RESISTANCE,";")))))),by="org"]
abricate.resfinder.tmp
lht.dt <- merge(lht.dt, abricate.resfinder.tmp, by="org")
# lht.dt
# cazymes
print("cazyme")
# substrate.hits: subfam hits with associated substrate
# substrate: total number of substrates
# substrate.unique: unique number of substrates
# Define cazyme classes
cazyme.classes=c("GH","GT","PL","CE","AA", "CBM")
print ('first')
dbcan.dt
# #colnames(lht.dt)
print("merge")
lht.dt <- merge(lht.dt, dbcan.dt[,list(dbcan_hits=.N, dbcan_consensus=sum(`#ofTools`==3), dbcan_signalp=sum(Signalp!="N")),by=org], by="org")
print((lht.dt))
print ('Second')
dbcan.dt
dbcan.tmp <- dcast(data=dbcan.dt[`#ofTools`==3,list(dbcan_class=str_extract(DIAMOND,paste0(cazyme.classes,collapse="|"))),by=org], formula=org~dbcan_class, fun=length)
dbcan.tmp
print("third")
lht.dt
print("dbcan.tmp")
print(dbcan.tmp)
print("colnames")
colnames(dbcan.tmp)
colnames(dbcan.tmp) <- ifelse(colnames(dbcan.tmp)=="org","org", paste0("dbcan_",colnames(dbcan.tmp)))
colnames(dbcan.tmp)
print(lht.dt$org)
print(dbcan.tmp$org)
lht.dt <- merge(lht.dt, dbcan.tmp, by="org")
lht.dt <- merge(lht.dt, dbcan.sub.dt[,list(dbcan_subfam.hits=.N, dbcan.substrates.hits=sum(Substrate!="-"), dbcan.substrates=sum(unlist(lapply(str_remove(unlist(str_split(Substrate,",")), "-"), function(x) str_length(x)>0))), dbcan.substate.unqiue=sum(unlist(lapply(unique(str_remove(trimws(unlist(str_split(Substrate,","))), "-")), function(x) str_length(x)>0)))),by=org], by="org")
head(lht.dt)
# antismash
print("antismash")
# antismash.dt[, org := as.character(org)]
# antismash.dt[, Type := as.character(Type)]
# antismash.dt[, Most.similar.known.cluster := as.character(Most.similar.known.cluster)]
# lht.dt <- merge(lht.dt, antismash.dt[,list(antismash.hits=.N, antismash.type=length(unique(Type)), antismash.known.cluster=sum(Most.similar.known.cluster!="")),by=org], by="org", all.x=T) # all.x=T because for some organisms no secondary metabolit cluster might found
# antismash.tmp <- dcast(data=antismash.dt[,list(type=unlist(str_split(Type, ","))), by=org], formula=org~type, fun=length)
# colnames(antismash.tmp) <- ifelse(colnames(antismash.tmp)=="org","org", paste0("antismash_",colnames(antismash.tmp)))
# lht.dt <- merge(lht.dt, antismash.tmp, by="org", all.x=T)
# lht.dt[is.na(lht.dt)] = 0 # set NA values from the merge of missing organisms to zero
# make sure org is character

antismash.dt[, org := as.character(org)]

if (nrow(antismash.dt) > 0) {

  lht.dt <- merge(
    lht.dt,
    antismash.dt[, .(
      antismash.hits = .N,
      antismash.type = length(unique(Type)),
      antismash.known.cluster = sum(Most.similar.known.cluster != "")
    ), by = org],
    by = "org",
    all.x = TRUE
  )

  antismash.tmp <- dcast(
    data = antismash.dt[
      , .(type = unlist(str_split(Type, ","), use.names = FALSE)),
      by = org
    ],
    formula = org ~ type,
    fun = length
  )

  colnames(antismash.tmp) <- ifelse(
    colnames(antismash.tmp) == "org",
    "org",
    paste0("antismash_", colnames(antismash.tmp))
  )

  lht.dt <- merge(lht.dt, antismash.tmp, by = "org", all.x = TRUE)

} else {
  # no hits → explicitly add zero columns
  lht.dt[, `:=`(
    antismash.hits = 0,
    antismash.type = 0,
    antismash.known.cluster = 0
  )]
}

# replace NAs introduced by merges
lht.dt[is.na(lht.dt)] <- 0


print("hua?")
lht.dt
# # gutsmash
print("gutsmash")
lht.dt <- merge(lht.dt, gutsmash.dt[,list(gutsmash.hits=.N, gutsmash.type=length(unique(Type)), gutsmash.known.cluster=sum(`Most similar known cluster`!="")),by=org], by="org", all.x=T)
lht.dt[is.na(lht.dt)] = 0 # set NA values from the merge of missing organisms to zero

# # # gapseq
print("gapseq")
metacyc.subsystems <- c("Activation-Inactivation-Interconversion","Bioluminescence","Biosynthesis","Degradation","Detoxification","Energy-Metabolism","Glycan-Pathways","Macromolecule-Modification","Metabolic-Clusters","Signaling-Pathways","Transport", "Enzyme-Test")

#db.meta <- rbind(fread("/home/garima/mlhp/data/gapseq/dat/meta_pwy.tbl"), fread("/home/garima/mlhp/data/gapseq/dat/custom_pwy.tbl"))
gapseq.pwy.dt <- merge(gapseq.pwy.dt, db.meta[,.(id,hierarchy)], by.x="ID", by.y="id")
gapseq.tmp <- dcast(data=gapseq.pwy.dt[Prediction==TRUE, list(subsystem=str_extract(hierarchy, paste0(metacyc.subsystems,collapse="|"))), by=org], formula=org~subsystem, fun=length)
colnames(gapseq.tmp) <- ifelse(colnames(gapseq.tmp)=="org","org", paste0("gapseq_meta.",colnames(gapseq.tmp)))
lht.dt <- merge(lht.dt, gapseq.tmp, by="org")
# print(head(lht.dt))
length(gapseq.models)
unique(names(gapseq.models))
# sapply(gapseq.models, class)
# lht.dt <- merge(lht.dt, data.table(org=names(gapseq.models), gapseq_growth=sapply(gapseq.models, function(x) sybil::optimizeProb(x)@lp_obj)), by="org")
lht.dt <- merge(lht.dt, data.table(org=names(gapseq.models), gapseq_growth=sapply(gapseq.models, function(x) cobrar::fba(x)@obj)), by="org")
lht.dt <- merge(lht.dt, gapseq.cs.dt[status==TRUE, list(gapseq_cs=.N), by=org], by="org")
lht.dt <- merge(lht.dt, gapseq.ferm.dt[status==TRUE, list(gapseq_ferm=.N), by=org], by="org")
print(head(lht.dt$gapseq_growth))
# eggnog
print("eggnog")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog.hits=.N, eggnog.kegg.tc=sum(unlist(lapply(str_remove(KEGG_TC, "-"), function(x) str_length(x)>0)))),by=org], by="org")
cog_category <- c <- c("J","A","K","L","B","D","Y","V","T","M","N","Z","W","U","O","C","G","E","F","H","I","P","Q","R","S")
eggnog.tmp <- dcast(data=eggnog.dt[,list(cog_category=unlist(str_extract_all(COG_category, paste0(cog_category, collapse="|")))),by=org], formula=org~cog_category, fun=length)
colnames(eggnog.tmp) <- ifelse(colnames(eggnog.tmp)=="org","org", paste0("eggnog_cog.cat.",colnames(eggnog.tmp)))
lht.dt <- merge(lht.dt, eggnog.tmp, by="org")
#eggnog.tmp2 <- dcast(data=eggnog.dt[,list(kegg_pwy=unlist(str_split(KEGG_Pathway,","))),by=org][kegg_pwy!="-"], formula=org~kegg_pwy, fun=length)
#colnames(eggnog.tmp2) <- ifelse(colnames(eggnog.tmp2)=="org","org", paste0("eggnog_kegg.pwy.",colnames(eggnog.tmp2)))
#lht.dt <- merge(lht.dt, eggnog.tmp2, by="org")
#eggnog.tmp3 <- dcast(data=eggnog.dt[,list(kegg_mod=unlist(str_split(KEGG_Module,","))),by=org][kegg_mod!="-"], formula=org~kegg_mod, fun=length)
#colnames(eggnog.tmp3) <- ifelse(colnames(eggnog.tmp3)=="org","org", paste0("eggnog_kegg.mod.",colnames(eggnog.tmp3)))
#lht.dt <- merge(lht.dt, eggnog.tmp3, by="org")
lht.dt
# platon
print("platon")
# platon.tmp <- platon.dt[,.(hits=.N, dbhits=sum(`# Plasmid Hits`), orfs=sum(`# ORFs`), amr=sum(`# AMRs`), replication=sum(`# Replication`),mobilization=sum(`# Mobilization`),conjugation=sum(`# Conjugation`), rrna=sum(`# rRNAs`), rds=mean(RDS)), by=org]
platon.tmp <- platon.dt[, .(
    hits = .N,                        # Count the number of rows for each group (i.e., how many entries per 'org')
    dbhits = sum(`# Plasmid Hits`),    # Sum the values in the column `# Plasmid Hits` for each group
    orfs = sum(`# ORFs`),              # Sum the values in the column `# ORFs` for each group
    amr = sum(`# AMRs`),               # Sum the values in the column `# AMRs` for each group
    replication = sum(`# Replication`),# Sum the values in the column `# Replication` for each group
    mobilization = sum(`# Mobilization`), # Sum the values in the column `# Mobilization` for each group
    conjugation = sum(`# Conjugation`),   # Sum the values in the column `# Conjugation` for each group
    rrna = sum(`# rRNAs`),             # Sum the values in the column `# rRNAs` for each group
    rds = mean(RDS)                    # Calculate the mean of the column `RDS` for each group
), by = org]
colnames(platon.tmp) <- ifelse(colnames(platon.tmp)=="org","org", paste0("platon_",colnames(platon.tmp)))
print("where is it starting from")
head(platon.tmp$org)
print("here comes lht.dt")
head(lht.dt$org)
platon.tmp[, org := as.character(org)]
lht.dt <- merge(lht.dt, platon.tmp, by="org", all.x=T)
print("really")
# Replace NAs (from empty platon or missing rows) with 0
lht.dt[is.na(lht.dt)] <- 0
# lht.dt[is.na(lht.dt)] = 0 # set NA values from the merge of missing organisms to zerok
head(lht.dt)
head(gapseq.med.dt$org)
merge(lht.dt, gapseq.med.dt[,list(gapseq_o2=sum(grepl("cpd00007",compounds))),by=org], by="org")
# aerobic lifestyle
lht.dt <- merge(lht.dt, gapseq.med.dt[,list(gapseq_o2=sum(grepl("cpd00007",compounds))),by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.aerobic.resp=sum(grepl("GO:0009060",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.anaerobic.resp=sum(grepl("GO:0009061",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.fermentation=sum(grepl("GO:0006113",GOs))), by=org], by="org")
lht.dt

print("stress")
# stress
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.stress.response=sum(grepl("GO:0006950",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.ros.response=sum(grepl("GO:0000302",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.o2.limit=sum(grepl("GO:0036293",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.o2.sensor=sum(grepl("GO:0019826",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.o2.response=sum(grepl("GO:0070482",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.ph.response=sum(grepl("GO:0009268",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.ph.regulation=sum(grepl("GO:0006885",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.ph.acidic=sum(grepl("GO:0010447",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.ph.alkaline=sum(grepl("GO:0010446",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.heat.response=sum(grepl("GO:0009408",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.phage.shock=sum(grepl("GO:0009271",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.oxidative.stress=sum(grepl("GO:0006979",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.osmotic.stress=sum(grepl("GO:0006970",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.dna.damage=sum(grepl("GO:0006974",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.water.deprivation=sum(grepl("GO:0009414",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.heavy.metal.response=sum(grepl("GO:0010038",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.efflux.export=sum(grepl("GO:0140352",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.efflux.multidrug.pump=sum(grepl("GO:0042910",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.membrane.protein.complex=sum(grepl("GO:0098797",GOs))), by=org], by="org")

# extra GO terms
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.thiamine.biosynthesis.pathways=sum(grepl("GO:0009228",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.osmoregulation=sum(grepl("GO:0009992",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.filamentous.growth=sum(grepl("GO:0030447",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.spore.dispersal=sum(grepl("GO:0075325",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.pigment.accumulation=sum(grepl("GO:0043476|GO:0042440",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.sporulation=sum(grepl("GO:0043934",GOs))), by=org], by="org")
# lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.anaerobic.respiration=sum(grepl("GO:0009061",GOs))), by=org], by="org")
# lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.aerobic.respiration=sum(grepl("GO:0009060",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.cellular.respiration=sum(grepl("GO:0045333",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.respiratory.gaseous.exchange=sum(grepl("GO:0007585",GOs))), by=org], by="org")

lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.amino.acid.biosynthetic=sum(grepl("GO:0008652",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.nucleotide.biosynthetic=sum(grepl("GO:0009165",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.fatty.acid.biosynthetic=sum(grepl("GO:0006633",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.translation=sum(grepl("GO:0006412",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.transport=sum(grepl("GO:0006810",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.methyl.accepting.chemotaxis.protein=sum(grepl("GO:0098561",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.mRNA.transcription=sum(grepl("GO:0009299",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.reactive.oxygen.species.biosynthetic=sum(grepl("GO:1903409",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.cell.growth=sum(grepl("GO:0016049",GOs))), by=org], by="org")

lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.membrane.lipid.metabolic=sum(grepl("GO:1905038",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.membrane.organisation=sum(grepl("GO:0061024",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.antioxidant.activity=sum(grepl("GO:0016209",GOs))), by=org], by="org")

# EPS
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.biofilm.attachment=sum(grepl("GO:0042711",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.exopolysaccharide.biosynthesis=sum(grepl("GO:0016209",GOs))), by=org], by="org")

# Dna repair
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.dna.repair=sum(grepl("GO:0006281",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.membrane.organisation=sum(grepl("GO:0061024",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.membrane.lipid.metabolic.regulation=sum(grepl("GO:1905038",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.carbepenim.biosynthesis=sum(grepl("GO:1901769",GOs))), by=org], by="org")

lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.free.radical.scavenger=sum(grepl("GO:0016209",GOs))), by=org], by="org")


# Stress - yas
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.response.to.desiccation=sum(grepl("GO:0009269",GOs))), by=org], by="org")
# lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.response.to.osmotic.stress=sum(grepl("GO:0006970",GOs))), by=org], by="org")  
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.trehalose.biosynthesis=sum(grepl("GO:0005992",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.glycine.betaine.biosynthesis=sum(grepl("GO:0031460",GOs))), by=org], by="org")  
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.double.strand.break.repair=sum(grepl("GO:0006302",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.protein.refolding=sum(grepl("GO:0042026",GOs))), by=org], by="org")  

lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.external.polysachharide.biosynthesis=sum(grepl("GO:0045226",GOs))), by=org], by="org")

lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.mismatch.repair=sum(grepl("GO:0006298",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.clavulanic.acid.biosynthesis=sum(grepl("GO:0033050",GOs))), by=org], by="org")  # enriched competative 
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.polyketide.metabolic=sum(grepl("GO:0030638",GOs))), by=org], by="org")  # competative dec
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.lipoprotein.transport=sum(grepl("GO:0042953",GOs))), by=org], by="org")  # competative enchanced 
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.biotin.metabolism=sum(grepl("GO:0006768",GOs))), by=org], by="org") # foraging 
                                                                              
# responsiveness
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.response.to.stimuli=sum(grepl("GO:0050896",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.cellular.response.to.chemical.stimuli=sum(grepl("GO:0070887",GOs))), by=org], by="org") 

# Extracellular enzymes (CSO stress)
# lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.proteolysis=sum(grepl("GO:0006508",GOs))), by=org], by="org") 
# low atp yielding metabolism
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.anaerobic.electron.transport=sum(grepl("GO:0019645",GOs))), by=org], by="org")  
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.oxidative.phosphorylation=sum(grepl("GO:0006119",GOs))), by=org], by="org")  


# extra GO terms till here!

# biofilm/attachement
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.biofilm.matrix=sum(grepl("GO:0062039",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.biofilm.formation=sum(grepl("GO:0042710",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.cell.adhesion=sum(grepl("GO:0007155",GOs))), by=org], by="org")

# # lipids
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.lipid.metaboism=sum(grepl("GO:0006629",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.sphingolipid.metabolism=sum(grepl("GO:0006665",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.sphingolipid.syn=sum(grepl("GO:0030148",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.sphingomyelin.metabolism=sum(grepl("GO:0006684",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.sphinganine.metabolism=sum(grepl("GO:0006668",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.glycosphingolipid.metabolism=sum(grepl("GO:0006687",GOs))), by=org], by="org")
gapseq.pwy.dt
lht.dt <- merge(lht.dt, dcast(data=gapseq.pwy.dt[grepl("Lipid-Biosynthesis", hierarchy),.(ID,Prediction), by=org], formula=org~ID, fun=sum, value.var="Prediction")[,list(gapseq_meta_lipid.syn=rowSums(.SD)),by=org,.SDcols=is.numeric], by="org")
lht.dt <- merge(lht.dt, dcast(data=gapseq.pwy.dt[grepl("Sphingolipid-Biosynthesis", hierarchy),.(ID,Prediction), by=org], formula=org~ID, fun=sum, value.var="Prediction")[,list(gapseq_meta_sphingolipid.syn=rowSums(.SD)),by=org,.SDcols=is.numeric], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.cholesterol.metabolism=sum(grepl("GO:0008203",GOs))), by=org], by="org")

# foraging
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.chemotaxis=sum(grepl("GO:0006935",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.motility=sum(grepl("GO:0048870",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_kegg.chemotaxis=sum(grepl("map02030",KEGG_Pathway))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_kegg.flagellar.assembly=sum(grepl("map02040",KEGG_Pathway))), by=org], by="org")

# signalling
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.quorum.sensing=sum(grepl("GO:0009372",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.signalling=sum(grepl("GO:0023052",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_kegg.two.component.system=sum(grepl("map02020",KEGG_Pathway))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_kegg.quorum.sensing=sum(grepl("map02024",KEGG_Pathway))), by=org], by="org")

# siderophores
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.siderophore.metabolism=sum(grepl("GO:0009237",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, gapseq.pwy.dt[Prediction==TRUE, list(gapseq_meta_siderophore.biosynthesis=length(grep("Siderophores-Biosynthesis", hierarchy))), by=org], by="org")

# spore/dormancy
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.dormancy=sum(grepl("GO:0022611",GOs))), by=org], by="org")

# antibiotics
lht.dt <- merge(lht.dt, gapseq.pwy.dt[Prediction==TRUE, list(gapseq_meta_antibiotic.biosynthesis=length(grep("Antibiotic-Biosynthesis", hierarchy))), by=org], by="org")
lht.dt <- merge(lht.dt, gapseq.pwy.dt[Prediction==TRUE, list(gapseq_meta_antibiotic.degradation=length(grep("ANTIBIOTIC-DEGRADATION", hierarchy))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_kegg.antibiotics.biosynthesis=sum(grepl("map00998",KEGG_Pathway))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.antibiotics.biosynthesis=sum(grepl("GO:0017000",GOs))), by=org], by="org")

# toxins
lht.dt <- merge(lht.dt, gapseq.pwy.dt[Prediction==TRUE, list(gapseq_meta_toxin.biosynthesis=length(grep("Toxin-Biosynthesis", hierarchy))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.toxin.biosynthesis=sum(grepl("GO:0009403",GOs))), by=org], by="org")

# bacteriocins
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.bacteriocin.response=sum(grepl("GO:0046678|GO:0071237",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.bacteriocin.transport=sum(grepl("GO:0043213",GOs))), by=org], by="org")

# sugar vs. acid catabolism
lht.dt <- merge(lht.dt, gapseq.pwy.dt[Prediction==TRUE, list(gapseq_meta_carboxylate.degradation=length(grep("CARBOXYLATES-DEG", hierarchy))), by=org], by="org")
lht.dt <- merge(lht.dt, gapseq.pwy.dt[Prediction==TRUE, list(gapseq_meta_fatty.acid.degradation=length(grep("Fatty-Acid-Degradation", hierarchy))), by=org], by="org")
lht.dt <- merge(lht.dt, gapseq.pwy.dt[Prediction==TRUE, list(gapseq_meta_sugar.degradation=length(grep("Sugars-And-Polysaccharides-Degradation", hierarchy))), by=org], by="org")

# proteolyse
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.proteolysis=sum(grepl("GO:0006508",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, gapseq.pwy.dt[Prediction==TRUE, list(gapseq_meta_h2s.biosynthesis=length(grep("Hydrogen-Sulfide-Biosynthesis", hierarchy))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.cresol.metabolism=sum(grepl("GO:0042212",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, gapseq.pwy.dt[Prediction==TRUE, list(gapseq_meta_aa.degradation=length(grep("Amino-Acid-Degradation", hierarchy))), by=org], by="org")

# nitrogen cycle
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.nitrogenous.reductase=sum(grepl("GO:0016661",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, gapseq.pwy.dt[Prediction==TRUE, list(gapseq_meta.nitrogen.metabolism=length(grep("NITROGEN-DEG", hierarchy))), by=org], by="org")
lht.dt <- merge(lht.dt, gapseq.pwy.dt[Prediction==TRUE, list(gapseq_meta.denitrification=length(grep("Denitrification", hierarchy))), by=org], by="org")

# secretion
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_kegg.secretion.systems=sum(grepl("map03070",KEGG_Pathway))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.protein.secretion=sum(grepl("GO:0009306",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.t1ss.secretion=sum(grepl("GO:0030253|GO:0030256",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.t2ss.secretion=sum(grepl("GO:0015628|GO:0015627",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.t3ss.secretion=sum(grepl("GO:0030254|GO:0030257",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.t4ss.secretion=sum(grepl("GO:0030255|GO:0043648",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.t5ss.secretion=sum(grepl("GO:0046819",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.t6ss.secretion=sum(grepl("GO:0033103|GO:0033104",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.extracellular.region=sum(grepl("GO:0005576",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.extracellular.space=sum(grepl("GO:0005615",GOs))), by=org], by="org")
print("try")
# print(lht.dt)
print(bakta.ext.dt)
lht.dt <- merge(lht.dt, bakta.ext.dt[,list(bakta_t1ss=sum(grepl("type i secretion",Product, ignore.case=T))), by=org])
lht.dt <- merge(lht.dt, bakta.ext.dt[,list(bakta_t2ss=sum(grepl("type ii secretion",Product, ignore.case=T))), by=org])
lht.dt <- merge(lht.dt, bakta.ext.dt[,list(bakta_t3ss=sum(grepl("type iii secretion",Product, ignore.case=T))), by=org])
lht.dt <- merge(lht.dt, bakta.ext.dt[,list(bakta_t4ss=sum(grepl("type iv secretion",Product, ignore.case=T))), by=org])
lht.dt <- merge(lht.dt, bakta.ext.dt[,list(bakta_t5ss=sum(grepl("type v secretion",Product, ignore.case=T))), by=org])
lht.dt <- merge(lht.dt, bakta.ext.dt[,list(bakta_t6ss=sum(grepl("type vi secretion",Product, ignore.case=T))), by=org])

# interaction
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.host.interactions=sum(grepl("GO:0051701",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.symbiotic.interactions=sum(grepl("GO:0044403",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.interactions=sum(grepl("GO:0044419",GOs))), by=org], by="org")

# necromass degradation, recycling
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.peptidoglycan.catabolism=sum(grepl("GO:0009253|GO:0009286",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, gapseq.pwy.dt[Prediction==TRUE, list(gapseq_meta.peptidoglycan.recycling=length(grep("Anhydromuropeptides-Recycling", hierarchy))), by=org], by="org")

# glycans, mucus etc
lht.dt <- merge(lht.dt, gapseq.pwy.dt[Prediction==TRUE, list(gapseq_meta.acgam.degradation=length(grep("N-Acetylglucosamine-Degradation", hierarchy))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.acgam.catabolism=sum(grepl("GO:0006046",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.murnac.catabolism=sum(grepl("GO:0097175|GO:0097174",GOs))), by=org], by="org") # N-acetyl-β-D-muramate
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.glucosamine.catabolism=sum(grepl("GO:1901072",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.chitin.catabolism=sum(grepl("GO:0006032",GOs))), by=org], by="org")

# regulation
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.sigma.factor=sum(grepl("GO:0016987",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.metabolic.regulation=sum(grepl("GO:0019222",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.gene.expression.gegulation=sum(grepl("GO:0010468",GOs))), by=org], by="org")

# polyamine
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_kegg.polyamine.biosynthesis=sum(grepl("M00133",KEGG_Module))), by=org], by="org")
lht.dt <- merge(lht.dt, gapseq.pwy.dt[Prediction==TRUE, list(gapseq_meta.amine.polyamine.biosynthesis=length(grep("Polyamine-Biosynthesis", hierarchy))), by=org], by="org")
lht.dt <- merge(lht.dt, gapseq.pwy.dt[Prediction==TRUE, list(gapseq_meta.amine.polyamine.degradation=length(grep("AMINE-DEG", hierarchy))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.polyamine.metabolism=sum(grepl("GO:0006595",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.polyamine.biosynthesis=sum(grepl("GO:0006596",GOs))), by=org], by="org")
lht.dt <- merge(lht.dt, eggnog.dt[,list(eggnog_go.polyamine.catabolism=sum(grepl("GO:0006598",GOs))), by=org], by="org")


write.csv(lht.dt, gzfile(paste0(input.dir, "/lht.csv.gz")), row.names=FALSE)
