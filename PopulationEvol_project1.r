library(VariantAnnotation)

# read vcf files
vcf_chr76M8M_YRI <- readVcf("chr7YRI_6000000_8000000.vcf", "hg19")
vcf_chr76M8M_CHB <- readVcf("chr7CHB_6000000_8000000.vcf", "hg19")

# row index to choose bialleclic snps
isSNV_biallelic_only_logic <- isSNV(vcf_chr76M8M_YRI,singleAltOnly=T)
isSNV_biallelic_only_logic_CHB <- isSNV(vcf_chr76M8M_CHB, singleAltOnly=T)


# subset of vcf with biallelic snps only
vcf_chr76M8M_YRI_snpsOnly <- vcf_chr76M8M_YRI[isSNV_biallelic_only_logic,]
vcf_chr76M8M_CHB_snpsOnly <- vcf_chr76M8M_CHB[isSNV_biallelic_only_logic_CHB,]


# extract geno "GT" information from the subset of vcf
vcf_chr76M8M_YRI_genoSimpleMatrix <- geno(vcf_chr76M8M_YRI_snpsOnly)[["GT"]]
vcf_chr76M8M_CHB_genoSimpleMatrix <- geno(vcf_chr76M8M_CHB_snpsOnly)[["GT"]]

detach("package:VariantAnnotation")
library(tidyverse)


# convert matrix to tibble
YRI_chr7_6M8Mgeno_tbl <- vcf_chr76M8M_YRI_genoSimpleMatrix %>% 
  as.data.frame(stringsAsFactors = F) %>% 
  rownames_to_column(var = "Id") %>% 
  as_tibble()
YRI_chr7_6M8Mgeno_tbl

CHB_chr7_6M8Mgeno_tbl <- vcf_chr76M8M_CHB_genoSimpleMatrix %>% 
  as.data.frame(stringsAsFactors = F) %>% 
  rownames_to_column(var = "Id") %>% 
  as_tibble()
CHB_chr7_6M8Mgeno_tbl


# convert to pivot-long tibble
YRI_chr7_6M8Mgeno_tbl_long <- YRI_chr7_6M8Mgeno_tbl %>% 
  pivot_longer(cols = -Id, names_to = "Samples", values_to = "GT")
YRI_chr7_6M8Mgeno_tbl_long 


CHB_chr7_6M8Mgeno_tbl_long <- CHB_chr7_6M8Mgeno_tbl %>% 
  pivot_longer(cols=-Id, names_to = "Samples", values_to = "GT")
CHB_chr7_6M8Mgeno_tbl_long


# add count columns to the tibble
YRI_chr7_6M8Mgeno_tbl_long_with_count <- YRI_chr7_6M8Mgeno_tbl_long %>% 
  mutate(homref_count = ifelse(GT == "0|0",1,0 )) %>%
  mutate(het_count = ifelse(GT == "0|1",1,0)) %>%
  mutate(homalt_count = ifelse(GT == "1|1",1,0)) %>%
  mutate(alt_count = homalt_count * 2 + het_count) %>%
  mutate(ref_count = homref_count * 2 + het_count )

YRI_chr7_6M8Mgeno_tbl_long_with_count


CHB_chr7_6M8Mgeno_tbl_long_with_count <- CHB_chr7_6M8Mgeno_tbl_long %>% 
  mutate(homref_count = ifelse(GT == "0|0",1,0 )) %>%
  mutate(het_count = ifelse(GT == "0|1",1,0)) %>%
  mutate(homalt_count = ifelse(GT == "1|1",1,0)) %>%
  mutate(alt_count = homalt_count * 2 + het_count) %>%
  mutate(ref_count = homref_count * 2 + het_count )  

YRI_chr7_6M8Mgeno_tbl_long_with_count

# create summary table 
YRI_chr7_6M8Mgeno_tbl_summary <- YRI_chr7_6M8Mgeno_tbl_long_with_count %>% 
  group_by(Id) %>% 
  summarise(alt_count_per_locus = sum(alt_count),
            ref_count_per_locus = sum(ref_count))

YRI_chr7_6M8Mgeno_tbl_summary


CHB_chr7_6M8Mgeno_tbl_summary <- CHB_chr7_6M8Mgeno_tbl_long_with_count %>% 
  group_by(Id) %>% 
  summarise(alt_count_per_locus = sum(alt_count),
            ref_count_per_locus = sum(ref_count))

CHB_chr7_6M8Mgeno_tbl_summary


# calculate minor allele frequencies  

YRI_chr7_6M8M_df <- YRI_chr7_6M8Mgeno_tbl_summary[,2:3] %>% data.frame()
YRI_chr7_6M8M_df$MAF <- apply(YRI_chr7_6M8M_df, 1,min)/apply(YRI_chr7_6M8M_df,1,sum) 
     
head(YRI_chr7_6M8M_df)  
table(round(YRI_chr7_6M8M_df$MAF,1))

CHB_chr7_6M8M_df <- CHB_chr7_6M8Mgeno_tbl_summary[,2:3] %>% data.frame()
CHB_chr7_6M8M_df$MAF <- apply(CHB_chr7_6M8M_df,1,min)/apply(CHB_chr7_6M8M_df,1,sum)

head(CHB_chr7_6M8M_df )
table(round(CHB_chr7_6M8M_df$MAF,1))

# histogram of SFS comparing YRI and CHB
MAF_df <- rbind(data.frame(population="YRI",MAF=YRI_chr7_6M8M_df$MAF),
                data.frame(population="CHB",MAF=CHB_chr7_6M8M_df$MAF))
head(MAF_df)

ggplot(MAF_df, aes(x=MAF, color=population))+
  geom_histogram()
  
  
