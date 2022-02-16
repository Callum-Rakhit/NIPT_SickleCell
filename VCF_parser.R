# Function to install packages
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

# Install required packages
GetPackages(c()) # currently all base R

# Get targets
targets <- read.delim(file = "/home/callum/snappy/snappy/tas/nipt/nipt-force-call-alleles.vcf", 
                      sep = "\t", header = F)
names(targets) <- as.character(unlist(targets[2,]))
targets <- targets[-c(1, 2),]

# Load VCF
VCF_Loader <- function(Input_VCF_Location){
  Input_VCF <- read.delim(file = Input_VCF_Location)
  colnames_vcf <- as.character(unlist(Input_VCF[128:137,]))
  Input_VCF <- Input_VCF[-c(1:137),]
  Output_VCF <- data.frame(matrix(ncol = 10, nrow = 0))
  colnames(Output_VCF) <- colnames_vcf
  x = 1
  y = 10
  for(i in 1:(length(Input_VCF)/10)){
    Output_VCF[nrow(Output_VCF) + 1,] <- as.character(unlist(t(Input_VCF[x:y])))
    x = x + 10
    y = y + 10
  }
  return(Output_VCF)
}

Final_Output <- targets[,1:5]

All_VCFs <- list.files(path = "/home/callum/Documents/NIPT_Sickle/NIPT_mutect_VCFs/", 
                       pattern = "*.vcf")
All_VCFs <- paste("/home/callum/Documents/NIPT_Sickle/NIPT_mutect_VCFs/", All_VCFs, sep = "")

Genotype_Extractor <- function(VCFfilepath){
  Output_VCF <- VCF_Loader(Input_VCF_Location = VCFfilepath)
  
  # Family 
  family_IDs <- read.delim(
    file = "/home/callum/Documents/NIPT_Sickle/Haplotyping_families.csv", 
    header = T, sep = ",")
  
  # Extract target genotypes from VCF
  
  # Get family information
  family_info <- family_IDs[grep(colnames(Output_VCF)[10], family_IDs$DNA.number), ]
  family_info <- paste(as.numeric(family_info$Family.ID), 
                       as.character(family_info$X), 
                       as.character(family_info$Genotype), 
                       as.character(family_info$DNA.number), 
                       sep = "_")
  
  # Get genotype information
  Output_Genotypes <- targets[,1:5]
  Output_Genotypes[,paste("Genotype", family_info, sep = ".")] <- NA
  Output_Genotypes[,paste("Allelic_Distribution", family_info, sep = ".")] <- NA
  
  # Collate together
  for(i in 1:length(Output_Genotypes$POS)){
    position <- as.character(unlist(targets[[2]][i]))
    temp_output <- Output_VCF[grep(position, Output_VCF$POS), ]
    if(dim(temp_output)!=0){
      Output_Genotypes[paste("Genotype", family_info, sep = ".")][i,] <- str_split(string = temp_output[,10], pattern = ":")[[1]][1]
      Output_Genotypes[paste("Allelic_Distribution", family_info, sep = ".")][i,] <- str_split(string = temp_output[,10], pattern = ":")[[1]][2]
    } else {
      Output_Genotypes[paste("Genotype", family_info, sep = ".")][i,] <- "Missing_genotype"
      Output_Genotypes[paste("Allelic_Distribution", family_info, sep = ".")][i,] <- "Missing_genotype"
    }
  }
  return(Output_Genotypes[,6:7])
}

temp <- lapply(All_VCFs, Genotype_Extractor)
for(i in 1:length(temp)){
  Final_Output <- cbind(Final_Output, temp[i])
}

Final_Output <- Final_Output[,order(colnames(Final_Output))]
write.table(x = Final_Output, 
            file = "/home/callum/Documents/NIPT_Sickle/Genotype_Output_v1.csv", 
            quote = F, col.names = T, row.names = F, sep = "\t")














