#install.package
#if (!requireNamespace("BiocManager", quietly = TRUE)) 
#  install.packages("BiocManager") 
#BiocManager::install("TCGAbiolinks")

#Load library 
library("TCGAbiolinks")
library(dplyr)
library(DT)
library(SummarizedExperiment)

setwd("/Users/kaurh8/Documents/GDC_TCGA_BiolinkS/")
a="TCGA"
b="LUSC"
# Create Query 
?GDCquery
query_TCGA_LUAD
paste0("query","_",a,"_",b)
query<- GDCquery(project = paste0(a,"-",b), 
                            data.category = "Transcriptome Profiling", 
                            data.type = "Gene Expression Quantification",  
                            # workflow.type = "HTSeq - FPKM", 
                            workflow.type = "STAR - Counts",
                            sample.type = "Primary Tumor" ) 



# look what type of information there in data 
LUAD_res = getResults(query_TCGA_LUAD) 



#query_TCGA_LUAD <- GDCquery(project = "TCGA-LUAD", data.category = "Transcriptome Profiling", data.type = "Gene Expression Quantification",  workflow.type = "HTSeq - FPKM", sample.type = c("Primary Tumor", "Solid Tissue Normal")) 
query_TCGA_LUAD <- GDCquery(project = "TCGA-LUAD", 
                            data.category = "Transcriptome Profiling", 
                            data.type = "Gene Expression Quantification",  
                           # workflow.type = "HTSeq - FPKM", 
                           workflow.type = "STAR - Counts",
                            sample.type = "Primary Tumor" ) 


# look what type of information there in data 

LUAD_res = getResults(query_TCGA_LUAD) 

#colnames(LUAD_res) 

head(LUAD_res$sample_type )

GDCdownload(query_TCGA_LUAD)


TCGA_LUAD_exp <- GDCprepare(query_TCGA_LUAD, 
                         save = TRUE, 
                         summarizedExperiment = TRUE, 
                         save.filename = "GDC_LUAD_Illumina_HiSeq.rda")
TCGA_LUAD_exp



?assays
TCGA_LUAD_Matrix <-  assay(TCGA_LUAD_exp) # only raw counts

TCGA_LUAD_Matrix_all <- assays(TCGA_LUAD_exp, withDimnames=TRUE) #for All types (raw counts.FPKM and TPM, unstranded(raw counts))
TCGA_LUAD_Matrix_all
#TCGA_LUAD_Matrix1 <- assays(TCGA_LUAD_exp,"stranded_first")
TCGA_LUAD_Matrix_raw_counts <- TCGA_LUAD_Matrix_all$unstranded
TCGA_LUAD_Matrix1_FPKM <- TCGA_LUAD_Matrix_all$fpkm_unstrand

TCGA_LUAD_Matrix_raw_counts1 <- as.data.frame(TCGA_LUAD_Matrix_raw_counts)
TCGA_LUAD_Matrix1_FPKM1 <- as.data.frame(TCGA_LUAD_Matrix1_FPKM )

write.table(cbind("ID"=rownames(TCGA_LUAD_Matrix_raw_counts1), TCGA_LUAD_Matrix_raw_counts1),file="TCGA_LUAD_Matrix_raw_counts.txt",sep="\t",quote=F, row.names=F)
write.table(cbind("ID"=rownames(TCGA_LUAD_Matrix1_FPKM1), TCGA_LUAD_Matrix1_FPKM1),file="TCGA_LUAD_Matrix_FPKM.txt",sep="\t",quote=F, row.names=F)

#gene annotation
TCGA_LUAD_gene_annotation <- as.data.frame(rowRanges(TCGA_LUAD_exp))
write.table(cbind("ID"=rownames(TCGA_LUAD_gene_annotation), TCGA_LUAD_gene_annotation),file="TCGA_LUAD_gene_annotation.txt",sep="\t",quote=F, row.names=F)


#clinincal data
clinical_TCGA_LUAD <- GDCquery_clinic(project = "TCGA-LUAD", type = "clinical")
clinical_TCGA_LUAD <- as.data.frame(clinical_TCGA_LUAD)
row.names(clinical_TCGA_LUAD) <- clinical_TCGA_LUAD$bcr_patient_barcode
write.table(cbind("ID"=rownames(clinical_TCGA_LUAD), clinical_TCGA_LUAD),file="clinical_TCGA_LUAD.txt",sep="\t",quote=F, row.names=F)


#more Clinical data
query_clin2 <- GDCquery(project = "TCGA-LUAD", 
                  data.category = "Clinical", 
                  file.type = "xml")

GDCdownload(query_clin2)

query_drug
?GDCprepare_clinic
clinical.drug <- GDCprepare_clinic(query_clin2, clinical.info = "drug")
clinical.drug <- as.data.frame(clinical.drug)
#row.names(clinical.drug) <- clinical.drug$bcr_patient_barcode
write.table(cbind("ID"=rownames(clinical.drug), clinical.drug),file="TCGA_LUAD_clinical.drug.txt",sep="\t",quote=F, row.names=F)

clinical.follow_up <- GDCprepare_clinic(query_clin2, clinical.info = "follow_up")
clinical.follow_up <- as.data.frame(clinical.follow_up)
#row.names(clinical.drug) <- clinical.drug$bcr_patient_barcode
write.table(cbind("ID"=rownames(clinical.follow_up), clinical.follow_up),file="TCGA_LUAD_clinical.follow_up.txt",sep="\t",quote=F, row.names=F)

clinical.NTE<- GDCprepare_clinic(query_clin2, clinical.info = "new_tumor_event")
clinical.NTE <- as.data.frame(clinical.NTE)
#row.names(clinical.drug) <- clinical.drug$bcr_patient_barcode
write.table(cbind("ID"=rownames(clinical.NTE), clinical.NTE),file="TCGA_LUAD_clinical.NTE.txt",sep="\t",quote=F, row.names=F)


clinical.pat <- GDCprepare_clinic(query_clin2, clinical.info = "patient")
clinical.pat <- as.data.frame(clinical.pat)
#row.names(clinical.drug) <- clinical.drug$bcr_patient_barcode
write.table(cbind("ID"=rownames(clinical.pat), clinical.pat),file="TCGA_LUAD_clinical.pat.txt",sep="\t",quote=F, row.names=F)


clinical.radiation <- GDCprepare_clinic(query_clin2, clinical.info = "radiation")
clinical.radiation <- as.data.frame(clinical.radiation)
#row.names(clinical.drug) <- clinical.drug$bcr_patient_barcode
write.table(cbind("ID"=rownames(clinical.radiation), clinical.radiation),file="TCGA_LUAD_clinical.radiation.txt",sep="\t",quote=F, row.names=F)

