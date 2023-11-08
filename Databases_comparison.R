# 8/11/2023
# Lucile Massenet-Regad


# Function for database comparison
comparison.lr.couple  <- function (db = db) {
  db.name.lr.couple = c()
  if (is.null(db$`Ligand 1`) | is.null(db$`Ligand 2`) | is.null(db$`Ligand 3`) |
      is.null(db$`Receptor 1`) | is.null(db$`Receptor 2`) | is.null(db$`Receptor 3`) | is.null(db$`Receptor 4`) | is.null(db$`Receptor 5`) ) {
    warning("Check database columns names : database should contains Ligand 1, Ligand 2, Ligand 3, Receptor 1, Receptor 2, Receptor 3, Receptor 4, and Receptor 5 column names")
  }
  for (mol in seq(1, dim(db)[1])) {
    #ligand name
    ligand_vector=c(db$`Ligand 1`[mol], db$`Ligand 2`[mol], db$`Ligand 3`[mol], db$`Ligand 4`[mol]) %>% na.omit() %>% as.vector()
    int=gtools::permutations(n=length(ligand_vector), r = length(ligand_vector), v =ligand_vector, repeats.allowed = F)
    ligand_name=apply(int,1,paste,collapse=" + ")
    
    # receptor name
    receptor_vector=c(db$`Receptor 1`[mol], db$`Receptor 2`[mol], db$`Receptor 3`[mol], db$`Receptor 4`[mol], db$`Receptor 5`[mol]) %>% na.omit() %>% as.vector()
    int2=gtools::permutations(n=length(receptor_vector), r = length(receptor_vector), v =receptor_vector, repeats.allowed = F)
    receptor_name=apply(int2,1,paste,collapse=" + ")
    
    # merge both
    names = c(do.call(paste0, expand.grid(ligand_name, " / ", receptor_name)), do.call(paste0, expand.grid(receptor_name, " / ", ligand_name)) )
    db.name.lr.couple = c(names, db.name.lr.couple)
  }
  db.name.lr.couple=unique(db.name.lr.couple)
  
  return(db.name.lr.couple)
}




###### ICELLNET DATABASE ####

## icellnet-db
DB=readxl::read_excel("ICELLNET/Databases/ICELLNETdb_V2.xlsx", guess_max = 3000)
DB$Database="ICELLNET"
DB$Reference=paste("PMID", DB$Reference)
DB=data.frame("Interaction_name" = NA, DB, check.names = F) # for interaction name to add in front
DB$Interaction_name=name.lr.couple(db=DB)[,1]



###### CellPhone DB ######
# https://github.com/ventolab/cellphonedb-data 

complex=read.csv(curl::curl(url="https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/complex_input.csv"), sep=",",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = "")
interaction=read.csv(curl::curl(url="https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/interaction_input.csv"), sep=",",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = "")
gene_info=read.csv(curl::curl(url="https://raw.githubusercontent.com/ventolab/cellphonedb-data/master/data/gene_input.csv"), sep=",",header = T, check.names=FALSE, stringsAsFactors = FALSE, na.strings = "")
gene_info=gene_info[!duplicated(gene_info$uniprot),] # remove duplicated rows for hgnc symbol and uniprot (only ensembl differs in these cases)

interaction= interaction %>% filter(annotation_strategy=="curated" & is_ppi =="True")
interaction$id_cp_interaction=paste0(interaction$partner_a, "_", interaction$partner_b)

db.new= as.data.frame(matrix(ncol=length(colnames(DB)) , nrow=length(interaction$id_cp_interaction)))
colnames(db.new)=colnames(DB)
db.new$Interaction_name=interaction$id_cp_interaction
db.new$Reference=interaction$source # source

# Transform cellphone db to be compatible with ICELLNET DB
q=c() # q <- store the vector of weird interaction that do not fit with ICELLNET format. 
for (mol in 1:length(db.new$Interaction_name)){
  #ligand
  if(interaction$partner_a[mol] %in% complex$complex_name){
    complex_info=complex[which(complex$complex_name==interaction$partner_a[mol]),]
    db.new$`Ligand 1`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_1)]
    if (!is.na(complex_info$uniprot_2)){
      db.new$`Ligand 2`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_2)]
    }
    if (!is.na(complex_info$uniprot_3)){
      db.new$`Ligand 3`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_3)]
    }
    if (!is.na(complex_info$uniprot_4)){
      db.new$`Ligand 4`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_4)]
    }
    if (!is.na(complex_info$uniprot_5)){
      q=c(q, mol)
    }
  }else{ db.new$`Ligand 1`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==interaction$partner_a[mol])] }
  
  #receptor
  if(interaction$partner_b[mol] %in% complex$complex_name){
    complex_info=complex[which(complex$complex_name==interaction$partner_b[mol]),]
    db.new$`Receptor 1`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_1)]
    if (!is.na(complex_info$uniprot_2)){
      db.new$`Receptor 2`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_2)]
    }
    if (!is.na(complex_info$uniprot_3)){
      db.new$`Receptor 3`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_3)]
    }
    if (!is.na(complex_info$uniprot_4)){
      db.new$`Receptor 4`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_4)]
    }
    if (!is.na(complex_info$uniprot_5)){
      db.new$`Receptor 5`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==complex_info$uniprot_5)]
    }
  }else{ db.new$`Receptor 1`[mol]=gene_info$hgnc_symbol[which(gene_info$uniprot==interaction$partner_b[mol])]}
}

db.new$Interaction_name=name.lr.couple(db=db.new)[,1]
#remove duplicate
db.new=db.new[!duplicated(db.new$Interaction_name),]




###### Cellchat DB ######

library(CellChat)
CellChatDB <- CellChatDB.human # see how add specific interactions
dim(CellChatDB$interaction)
head(CellChatDB$complex)
head(CellChatDB$cofactor) #not taken into account in ICELLNET
head(CellChatDB$geneInfo)

interaction=CellChatDB$interaction
complex=CellChatDB$complex
complex <- replace(complex, complex=="", NA)
head(complex)
complex$complex_name=rownames(complex)

db.new2= as.data.frame(matrix(ncol=length(colnames(DB)) , nrow=length(CellChatDB$interaction$interaction_name)))
colnames(db.new2)=colnames(DB)
db.new2$Interaction_name=CellChatDB$interaction$interaction_name
db.new2$Reference=CellChatDB$interaction$evidence #source


# Transform into format compatible with ICELLNET DB
for (mol in 1:length(db.new2$Interaction_name)){
  #mol=1
  #ligand
  if(interaction$ligand[mol] %in% complex$complex_name){
    db.new2$`Ligand 1`[mol]=complex$subunit_1[which(complex$complex_name==interaction$ligand[mol])]
    if(!is.na(complex$subunit_2[which(complex$complex_name==interaction$ligand[mol])])){
      db.new2$`Ligand 2`[mol]=complex$subunit_2[which(complex$complex_name==interaction$ligand[mol])]
      if (!is.na(complex$subunit_3[which(complex$complex_name==interaction$ligand[mol])])){
        db.new2$`Ligand 3`[mol]=complex$subunit_3[which(complex$complex_name==interaction$ligand[mol])]
        if (!is.na(complex$subunit_4[which(complex$complex_name==interaction$ligand[mol])])){
          db.new2$`Ligand 4`[mol]=complex$subunit_4[which(complex$complex_name==interaction$ligand[mol])]
        }
      }
    }
  }else{ db.new2$`Ligand 1`[mol]=interaction$ligand[mol] }
  
  #receptor
  if(interaction$receptor[mol] %in% complex$complex_name){
    db.new2$`Receptor 1`[mol]=complex$subunit_1[which(complex$complex_name==interaction$receptor[mol])]
    if(!is.na(complex$subunit_2[which(complex$complex_name==interaction$receptor[mol])])){
      db.new2$`Receptor 2`[mol]=complex$subunit_2[which(complex$complex_name==interaction$receptor[mol])]
      if (!is.na(complex$subunit_3[which(complex$complex_name==interaction$receptor[mol])])){
        db.new2$`Receptor 3`[mol]=complex$subunit_3[which(complex$complex_name==interaction$receptor[mol])]
        if (!is.na(complex$subunit_4[which(complex$complex_name==interaction$receptor[mol])])){
          db.new2$`Receptor 4`[mol]=complex$subunit_4[which(complex$complex_name==interaction$receptor[mol])]
        }
      }
    }
  }else{ db.new2$`Receptor 1`[mol]=interaction$receptor[mol] }
}

db.new2$Interaction_name=name.lr.couple(db=db.new2)[,1]