#install.package
if (!requireNamespace("BiocManager", quietly = TRUE)) 
  install.packages("BiocManager") 
BiocManager::install("TCGAbiolinks")

#Load library 
library("TCGAbiolinks")
library(dplyr)
library(DT)
library(SummarizedExperiment)


#setwd("/Users/kaurh8/Documents/GDC_TCGA_Biolinks/")


path="/Users/kaurh8/Documents/GDC_TCGA_Biolinks/"
cancer_list <- read.table("cancer_list",header =T, sep = "\t",   row.names = 1, check.names = FALSE)
cancer_list
list1 <- as.character(row.names(cancer_list))
list1

for (c in 1:length(list1))
  
{
  
  cancer <-  paste0(list1[c])
  print(cancer)
  paste0(path,list1[c], "/")
  
  dir.create(paste0(path,list1[c], "/"))
  setwd(paste0(path,list1[c], "/"))
  getwd()
  
  path_c <- paste0(path,list1[c], "/")
  #print("cancer folder", path_c)
  print(path_c)
  
  setwd(path_c)
  
  getwd()

  tryCatch({

    
query<- GDCquery(project = cancer, 
                 data.category = "Transcriptome Profiling", 
                 data.type = "Gene Expression Quantification",  
                 # workflow.type = "HTSeq - FPKM", 
                 workflow.type = "STAR - Counts",
                 sample.type = "Primary Tumor" ) 



# look what type of information there in data 
cancer_res = getResults(query) 


GDCdownload(query)


#TCGA_exp <- GDCprepare(query,  save = TRUE, 
#                            summarizedExperiment = TRUE, 
#                            save.filename = "GDC_LUAD_Illumina_HiSeq.rda")
#                             save.filename = paste0(i,"_Illumina_HiSeq.rda"))

TCGA_exp <- GDCprepare(query,  save = TRUE,  
                       summarizedExperiment = TRUE, 
                       save.filename = paste0(cancer,"_Illumina_HiSeq.rda")) 
                       

                      
library(SummarizedExperiment)


#?assays
#TCGA_LUAD_Matrix <-  assay(TTCGA_exp) # only raw counts

TCGA_Matrix_all <- assays(TCGA_exp, withDimnames=TRUE) #for All types (raw counts.FPKM and TPM, unstranded(raw counts))
TCGA_Matrix_all
#TCGA_LUAD_Matrix1 <- assays(TCGA_LUAD_exp,"stranded_first")
TCGA_Matrix_raw_counts <- TCGA_Matrix_all$unstranded
TCGA_Matrix1_FPKM <- TCGA_Matrix_all$fpkm_unstrand

TCGA_Matrix_raw_counts1 <- as.data.frame(TCGA_Matrix_raw_counts)
TCGA_Matrix1_FPKM1 <- as.data.frame(TCGA_Matrix1_FPKM )

write.table(cbind("ID"=rownames(TCGA_Matrix_raw_counts1), TCGA_Matrix_raw_counts1),file=paste0(cancer,"_Matrix_raw_counts.txt"),sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(TCGA_Matrix1_FPKM1), TCGA_Matrix1_FPKM1),file=paste0(cancer,"_Matrix_FPKM.txt"),sep="\t",quote=F, row.names=F)

#gene annotation
TCGA_gene_annotation <- as.data.frame(rowRanges(TCGA_exp))
write.table(cbind("ID"=rownames(TCGA_gene_annotation), TCGA_gene_annotation),file=paste0(cancer,"_gene_annotation.txt"),sep="\t",quote=F, row.names=F)


#clinincal data
clinical_TCGA <- GDCquery_clinic(project = cancer, type = "clinical")
clinical_TCGA <- as.data.frame(clinical_TCGA)
row.names(clinical_TCGA) <- clinical_TCGA$bcr_patient_barcode
write.table(cbind("ID"=rownames(clinical_TCGA), clinical_TCGA),file=paste0(cancer, "_clinical.txt"),sep="\t",quote=F, row.names=F)


#more Clinical data
query_clin2 <- GDCquery(project = cancer, 
                        data.category = "Clinical", 
                        file.type = "xml")

GDCdownload(query_clin2)


clinical.drug <- GDCprepare_clinic(query_clin2, clinical.info = "drug")
clinical.drug <- as.data.frame(clinical.drug)
#row.names(clinical.drug) <- clinical.drug$bcr_patient_barcode
write.table(cbind("ID"=rownames(clinical.drug), clinical.drug),file=paste0(cancer,"_clinical_drug.txt"),sep="\t",quote=F, row.names=F)

clinical.follow_up <- GDCprepare_clinic(query_clin2, clinical.info = "follow_up")
clinical.follow_up <- as.data.frame(clinical.follow_up)
#row.names(clinical.drug) <- clinical.drug$bcr_patient_barcode
write.table(cbind("ID"=rownames(clinical.follow_up), clinical.follow_up),file=paste0(cancer,"_clinical_followup.txt"),sep="\t",quote=F, row.names=F)

clinical.NTE<- GDCprepare_clinic(query_clin2, clinical.info = "new_tumor_event")
clinical.NTE <- as.data.frame(clinical.NTE)
#row.names(clinical.drug) <- clinical.drug$bcr_patient_barcode
write.table(cbind("ID"=rownames(clinical.NTE), clinical.NTE),file=paste0(cancer,"_clinical_NTE.txt"), sep="\t",quote=F, row.names=F)


clinical.pat <- GDCprepare_clinic(query_clin2, clinical.info = "patient")
clinical.pat <- as.data.frame(clinical.pat)
#row.names(clinical.drug) <- clinical.drug$bcr_patient_barcode
write.table(cbind("ID"=rownames(clinical.pat), clinical.pat),file=paste0(cancer,"_clinical_patients.txt") ,sep="\t",quote=F, row.names=F)


clinical.radiation <- GDCprepare_clinic(query_clin2, clinical.info = "radiation")
clinical.radiation <- as.data.frame(clinical.radiation)
#row.names(clinical.drug) <- clinical.drug$bcr_patient_barcode
write.table(cbind("ID"=rownames(clinical.radiation), clinical.radiation),file=paste0(cancer,"_clinical_radiation_trt.txt"),sep="\t",quote=F, row.names=F)

  }, error=function(e){cat("ERROR :",conditionMessage(e), "\n")})
}


