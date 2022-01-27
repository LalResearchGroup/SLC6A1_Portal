################################################
################# SLC6A1 PORTAL ################
################################################
################## DATA PREPARATION ############

library(readr)
library(tidyverse)


# Patient variants table
## read in old table for comparison
pat.dat.old <- read_delim("data/Patient_variants_SLC6A1_v5.txt", 
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