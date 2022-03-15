# Function to install packages
GetPackages <- function(required.packages) {
  packages.not.installed <- required.packages[!(required.packages %in% installed.packages()[, "Package"])]
  if(length(packages.not.installed)){install.packages(packages.not.installed, dependencies = T)}
  suppressMessages(lapply(required.packages, require, character.only = T))
}

# Install required packages
GetPackages(c("devtools", "stringr", "openxlsx", "XLConnect"))

# Can't have xlsx and XLConnect at the same time
detach(name = "xlsx", unload = T)
install_version("XLConnectJars", version = "0.2-12", repos = "http://cran.us.r-project.org") 
install_version("XLConnect", version = "0.2-12", repos = "http://cran.us.r-project.org")

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
All_VCFs <- list.files(path = "/home/callum/Documents/NIPT_Sickle/NIPT_mutect_VCFs/NIPT_mutect_VCFs_v2/", 
                       pattern = "*.vcf")
All_VCFs <- paste("/home/callum/Documents/NIPT_Sickle/NIPT_mutect_VCFs/NIPT_mutect_VCFs_v2/", All_VCFs, sep = "")

Genotype_Extractor <- function(VCFfilepath){
  Output_VCF <- VCF_Loader(Input_VCF_Location = VCFfilepath)
  # Output_VCF <- VCF_Loader(Input_VCF_Location = All_VCFs[1])
  
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
  Output_Genotypes[,paste(family_info, "a_Allelic_Distribution", sep = ".")] <- NA
  Output_Genotypes[,paste(family_info, "b_VCF_Reference_Genotype", sep = ".")] <- NA
  Output_Genotypes[,paste(family_info, "c_VCF_Alternative_Genotype", sep = ".")] <- NA
  Output_Genotypes[,paste(family_info, "d_Simple_Genotype", sep = ".")] <- NA
  Output_Genotypes[,paste(family_info, "e_Allele_Frequency", sep = ".")] <- NA
  Output_Genotypes[,paste(family_info, "f_Allele_1", sep = ".")] <- NA
  Output_Genotypes[,paste(family_info, "g_Allele_2", sep = ".")] <- NA
  
  # Collate together
  options(scipen = 999)
  
  for(i in 1:length(Output_Genotypes$POS)){
    position <- as.character(unlist(targets[[2]][i]))
    temp_output <- Output_VCF[grep(position, Output_VCF$POS), ]
    if(dim(temp_output)!=0){
      Output_Genotypes[paste(family_info, "a_Allelic_Distribution", sep = ".")][i,] <- str_split(string = temp_output[,10], pattern = ":")[[1]][2]
      Output_Genotypes[paste(family_info, "b_VCF_Reference_Genotype", sep = ".")][i,] <- temp_output[4]
      Output_Genotypes[paste(family_info, "c_VCF_Alternative_Genotype", sep = ".")][i,] <- temp_output[5]
      
      genotype = str_split(string = temp_output[,10], pattern = ":")[[1]][2]
      genotype_split = str_split(string = genotype, pattern = ",")
      all_reads = sum(as.numeric(genotype_split[[1]][1:length(genotype_split[[1]])]))
      reference = as.numeric(genotype_split[[1]][1])
      all_alternate = sum(as.numeric(genotype_split[[1]][2:length(genotype_split[[1]])]))
      
      if(all_reads < 10) {
        Output_Genotypes[paste(family_info, "d_Simple_Genotype", sep = ".")][i,] <- "Too_few_reads"
      } else if(reference == 0) {
        Output_Genotypes[paste(family_info, "d_Simple_Genotype", sep = ".")][i,] <- "Homozygous_Alternative(s)"
      } else if(all_alternate == 0) {
        Output_Genotypes[paste(family_info, "d_Simple_Genotype", sep = ".")][i,] <- "Homozygous_Reference"
      } else if(reference/all_alternate > 1.8) {
        Output_Genotypes[paste(family_info, "d_Simple_Genotype", sep = ".")][i,] <- "Homozygous_Reference"
      } else if(reference/all_alternate > 0.2) {
        Output_Genotypes[paste(family_info, "d_Simple_Genotype", sep = ".")][i,] <- "Hetrozygous"
      } else {
        Output_Genotypes[paste(family_info, "d_Simple_Genotype", sep = ".")][i,] <- "Homozygous_Alternative(s)"
      }
      
      Output_Genotypes[paste(family_info, "e_Allele_Frequency", sep = ".")][i,] <- str_split(string = temp_output[,10], pattern = ":")[[1]][3]
      Output_Genotypes[paste(family_info, "f_Allele_1", sep = ".")][i,] <- temp_output[4]
      
      if(all_reads > 10) {
        read_counts_temp <- str_split(string = str_split(string = temp_output[,10], pattern = ":")[[1]][2], pattern = ",")
        highest_value <- match(max(as.numeric(read_counts_temp[[1]][2:length(read_counts_temp[[1]])])), read_counts_temp[[1]][2:length(read_counts_temp[[1]])])
        Output_Genotypes[paste(family_info, "g_Allele_2", sep = ".")][i,] <- str_split(string = temp_output[5], pattern = ",")[[1]][highest_value]
      }
      
    } else {
      Output_Genotypes[paste(family_info, "a_Allelic_Distribution", sep = ".")][i,] <- "Missing_genotype"
      Output_Genotypes[paste(family_info, "b_VCF_Reference_Genotype", sep = ".")][i,] <- "Missing_genotype"
      Output_Genotypes[paste(family_info, "c_VCF_Alternative_Genotype", sep = ".")][i,] <- "Missing_genotype"
      Output_Genotypes[paste(family_info, "d_Simple_Genotype", sep = ".")][i,] <- "Missing_genotype"
      Output_Genotypes[paste(family_info, "e_Allele_Frequency", sep = ".")][i,] <- "Missing_genotype" 
      Output_Genotypes[paste(family_info, "f_Allele_1", sep = ".")][i,] <- "Missing_genotype"
      Output_Genotypes[paste(family_info, "g_Allele_2", sep = ".")][i,] <- "Missing_genotype"
    }
  }
  return(Output_Genotypes[,6:12])
}

temp <- lapply(All_VCFs, Genotype_Extractor)

for(i in 1:length(temp)){
  Final_Output <- cbind(Final_Output, temp[i])
}

# Remove missing genotypes
# Final_Output <- Final_Output[!grepl("Missing_genotype", Final_Output$`19_Mother_AC_14-1981.a_Allelic_Distribution`),]
Final_Output <- Final_Output[,order(colnames(Final_Output))]

# Arrange df variables by position, variables must be a named vector, e.g. c("var.name" = sensible.name)
arrange.vars <- function(data, vars){
  # Stop if not a data.frame (but should work for matrices as well)
  stopifnot(is.data.frame(data))
  
  # Sort out inputs
  data.nms <- names(data)
  var.nr <- length(data.nms)
  var.nms <- names(vars)
  var.pos <- vars
  
  # Sanity checks
  stopifnot( !any(duplicated(var.nms)), 
             !any(duplicated(var.pos)) )
  stopifnot( is.character(var.nms), 
             is.numeric(var.pos) )
  stopifnot( all(var.nms %in% data.nms) )
  stopifnot( all(var.pos > 0), 
             all(var.pos <= var.nr) )
  
  # Prepare output
  out.vec <- character(var.nr)
  out.vec[var.pos] <- var.nms
  out.vec[-var.pos] <- data.nms[ !(data.nms %in% var.nms) ]
  stopifnot( length(out.vec)==var.nr )
  
  # Rearrange variables by position
  data <- data[ , out.vec]
  return(data)
}

Final_Output <- arrange.vars(Final_Output, c("#CHROM" = 1, "ID" = 2, "POS" = 3, "REF" = 4, "ALT" = 5))

# Write as Excel with each family in a separate worksheet
file <- '/home/callum/Documents/NIPT_Sickle/Genotype_Output_v6.xlsx'

wb <- XLConnect::loadWorkbook(filename = file, create = T)
for(i in 1:20){
  sheet <- XLConnect::createSheet(object = wb, name = paste("Family_", i, sep = ""))
  family_sheet <- cbind(Final_Output[,1:5], Final_Output[grepl(paste("^", i, "_", sep = ""), names(Final_Output))])
  appendWorksheet(object = wb, data = family_sheet, sheet = paste("Family_", i, sep = ""), 
                  header = T, rownames = F)
}
XLConnect::saveWorkbook(object = wb, file = file)

write.table(x = Final_Output, 
            file = "/home/callum/Documents/NIPT_Sickle/Genotype_Output_v5.csv", 
            quote = T, col.names = T, row.names = F, sep = "\t")

#### Test section #####

Output_VCF <- VCF_Loader(Input_VCF_Location = All_VCFs[1])

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
Output_Genotypes[,paste(family_info, "a_Allelic_Distribution", sep = ".")] <- NA
Output_Genotypes[,paste(family_info, "b_VCF_Reference_Genotype", sep = ".")] <- NA
Output_Genotypes[,paste(family_info, "c_VCF_Alternative_Genotype", sep = ".")] <- NA
Output_Genotypes[,paste(family_info, "d_Simple_Genotype", sep = ".")] <- NA
Output_Genotypes[,paste(family_info, "e_Allele_Frequency", sep = ".")] <- NA
Output_Genotypes[,paste(family_info, "f_Allele_1", sep = ".")] <- NA
Output_Genotypes[,paste(family_info, "g_Allele_2", sep = ".")] <- NA

# Collate together
options(scipen = 999)

for(i in 1:length(Output_Genotypes$POS)){
  position <- as.character(unlist(targets[[2]][i])) # Needs [[2]][i]
  temp_output <- Output_VCF[grep("4904467", Output_VCF$POS), ]
  if(dim(temp_output)!=0){
    Output_Genotypes[paste(family_info, "a_Allelic_Distribution", sep = ".")][i,] <- str_split(string = temp_output[,10], pattern = ":")[[1]][2]
    Output_Genotypes[paste(family_info, "b_VCF_Reference_Genotype", sep = ".")][i,] <- temp_output[4]
    Output_Genotypes[paste(family_info, "c_VCF_Alternative_Genotype", sep = ".")][i,] <- temp_output[5]
    
    genotype = str_split(string = temp_output[,10], pattern = ":")[[1]][2]
    genotype_split = str_split(string = genotype, pattern = ",")
    all_reads = sum(as.numeric(genotype_split[[1]][1:length(genotype_split[[1]])]))
    reference = as.numeric(genotype_split[[1]][1])
    all_alternate = sum(as.numeric(genotype_split[[1]][2:length(genotype_split[[1]])]))
    
    if(all_reads < 10) {
      Output_Genotypes[paste(family_info, "d_Simple_Genotype", sep = ".")][i,] <- "Too_few_reads"
    } else if(reference == 0) {
      Output_Genotypes[paste(family_info, "d_Simple_Genotype", sep = ".")][i,] <- "Homozygous_Alternative(s)"
    } else if(all_alternate == 0) {
      Output_Genotypes[paste(family_info, "d_Simple_Genotype", sep = ".")][i,] <- "Homozygous_Reference"
    } else if(reference/all_alternate > 1.8) {
      Output_Genotypes[paste(family_info, "d_Simple_Genotype", sep = ".")][i,] <- "Homozygous_Reference"
    } else if(reference/all_alternate > 0.2) {
      Output_Genotypes[paste(family_info, "d_Simple_Genotype", sep = ".")][i,] <- "Hetrozygous"
    } else {
      Output_Genotypes[paste(family_info, "d_Simple_Genotype", sep = ".")][i,] <- "Homozygous_Alternative(s)"
    }
    
    Output_Genotypes[paste(family_info, "e_Allele_Frequency", sep = ".")][i,] <- str_split(string = temp_output[,10], pattern = ":")[[1]][3]
    Output_Genotypes[paste(family_info, "f_Allele_1", sep = ".")][i,] <- temp_output[4]
    
    if(all_reads > 10) {
      read_counts_temp <- str_split(string = str_split(string = temp_output[,10], pattern = ":")[[1]][2], pattern = ",")
      highest_value <- match(max(as.numeric(read_counts_temp[[1]][2:length(read_counts_temp[[1]])])), read_counts_temp[[1]][2:length(read_counts_temp[[1]])])
      Output_Genotypes[paste(family_info, "g_Allele_2", sep = ".")][i,] <- str_split(string = temp_output[5], pattern = ",")[[1]][highest_value]
    }
    
  } else {
    Output_Genotypes[paste(family_info, "a_Allelic_Distribution", sep = ".")][i,] <- "Missing_genotype"
    Output_Genotypes[paste(family_info, "b_VCF_Reference_Genotype", sep = ".")][i,] <- "Missing_genotype"
    Output_Genotypes[paste(family_info, "c_VCF_Alternative_Genotype", sep = ".")][i,] <- "Missing_genotype"
    Output_Genotypes[paste(family_info, "d_Simple_Genotype", sep = ".")][i,] <- "Missing_genotype"
    Output_Genotypes[paste(family_info, "e_Allele_Frequency", sep = ".")][i,] <- "Missing_genotype" 
    Output_Genotypes[paste(family_info, "f_Allele_1", sep = ".")][i,] <- "Missing_genotype"
    Output_Genotypes[paste(family_info, "g_Allele_2", sep = ".")][i,] <- "Missing_genotype"
  }
}
