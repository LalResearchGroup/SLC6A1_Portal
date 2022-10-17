################################################
################# SLC6A1 PORTAL ################
################################################
################## DATA PREPARATION ############

# Packages

library(readr)
library(tidyverse)
library(here)
library(readxl)

# Pat variant table
Pat_var.df <- read_delim("prep_data/Patient_variants_SLC6A1_v6_fixed.txt", "\t") %>% 
  select(-X27)
# Updated table as of October 10, 2022
Pat_var.df <- read_delim("prep_data/Patient_variants_SLC6A1_v7_edit.txt", "\t")

# output new file into 'data' table for use in ui/server
write_delim(Pat_var.df, "data/Patient_variants_SLC6A1_v7.txt", "\t") # latest file as of October 17, 2022 is v8 (patient removed from cohort due to lacking consent)

# prep new functional data from biomarin
df_biomarin_results_source <- read_excel(here("prep_data", "SLC6A1.supplement.tables.v03.xlsx"), 
           sheet = "4-Results") %>% 
  mutate(Gene = "SLC6A1") %>% 
  mutate(AA_pos = str_extract(HGVSp,"[0-9]+") %>% as.numeric()) %>% 
  mutate(AA_ref = str_extract(HGVSp, "^\\D+")) %>% 
  mutate(AA_ref = gsub("^.*?\\.","", AA_ref)) %>% 
  mutate(AA_alt = ifelse(`Variant impact` == "missense", sub(".*[0-9]", "", HGVSp), "NA")) %>% 
  mutate(variant_impact = case_when(AA_alt == "*" ~ "stop gained",
                                    AA_alt != "NA" ~ "missense",
                                    `Variant impact` == "frameshift" ~ "frameshift",
                                    `Variant impact` == "stop gained" ~ "stop gained",
                                    `Variant impact` == "inframe indel" ~ "inframe indel",
                                    `Variant impact` == "synonymous" ~ "synonymous")) %>% 
  mutate(uptake = Avg_PercentWT/100)
  #right_join(Pat_var.df, by = c("Pos" = "cDNA_pos", "Ref" = "cDNA_ref", "Alt" = "cDNA_alt")) %>% 
  #mutate(cDNA_pos = str_extract(Transcript,"[0-9]+") %>% as.numeric())

  
# output new functional data file for use in ui/server
write_csv(df_biomarin_results_source, "data/Functional_data_biomarin.csv")

# read in new domain data
df_domain <- read_delim(here("prep_data", "Domain_cornelius.txt"), delim = "\t") %>% 
  rename(Domain_old = Domain, Domain_color_old = Domain_color)
# output updated file to data
write_delim(df_domain, here("data", "Domain_cornelius.txt"), delim = "\t")


















################################################
################################################
############ OLD pat var table changes #########
##### check archived files for backtrace #######
################################################
################################################
# Patient variants table
## read in old table for comparison
pat.dat.old <- read_delim("data_archive/Patient_variants_SLC6A1_v5.txt", 
                          "\t", escape_double = FALSE, col_types = cols(Chr = col_character(), 
                                                                        X27 = col_character(), `Cognitive Impairment` = col_character()), 
                          trim_ws = TRUE) %>% 
  mutate(Genomic_pos = as.numeric(Genomic_pos), `Age at Inclusion (years)` = as.numeric(`Age at Inclusion (years)`))
# select(Transcript:cDNA_alt, AA_pos:AA_alt)
## save as csv
write_csv(pat.dat.old, "prep_data/Patient_variants_SLC6A1.csv")

# ## read in old mod table
# Patient_data.df.mod <- read_delim("data/Patient_variants_SLC6A1_mod.txt", delim = "\t")

## read in new v4 table
# pat.dat.v4 <- read_csv("prep_data/Patient_variants_SLC6A1_v4.csv") # take from prep_data folder
pat.dat.v4 <- read_excel("prep_data/Patient_variants_SLC6A1_v4new.xlsx", 
                                                          col_types = c("text", "text", "text", 
                                                                        "numeric", "numeric", "text", "text", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "numeric", "numeric", "numeric", 
                                                                        "text", "text", "text", "numeric", 
                                                                        "text", "text", "text", "text", "text", 
                                                                        "text", "text", "numeric", "text", 
                                                                        "text", "text", "numeric", "text", 
                                                                        "text", "text", "text", "text", "text", 
                                                                        "text", "text", "text", "text", "text", 
                                                                        "text", "text", "text", "text", "text", 
                                                                        "text", "text", "numeric", "text", 
                                                                        "text", "numeric"))


## prep for column swap with Tobias' updated table 'Patient_variants_SLC6A1.txt'
### join new columns from Tobias' table (v5) with new v4 table
pat.dat.v4 <- pat.dat.v4 %>% 
  select(-AA_pos, -AA_alt, -AA_ref) %>% 
  left_join(pat.dat.old) # this join removes entries in column X27, idk why!???

## save new v4 table as .txt tab del
write_delim(pat.dat.v4, "data/Patient_variants_SLC6A1_v4.txt", # write to data folder
            delim = "\t")
# ## save as csv
# write_csv(pat.dat.v4, "prep_data/Patient_variants_SLC6A1_v4test.csv")


#########################################
########### TRY FINDING PLOT ############

Patient_data.df