library(tidyverse)
library(Rfast)


#Calculate the variant burden inside a 3D bubble.

# Current script considers the single subunits only (GRIN1/ GRIN2A) will be explored in on subunit level rather than the protein complex. The gene for which the 3D variant burden will be calculated has to be chosen in advance. 
# If desired protein complex burden can be added. 


#radius for 3D bubble 
radius_selected = 12 # radius in Angstrom

#Required functions

min_X <- function(vec,threshold){
  
  return(order(vec)[1:threshold] %>% list())
  
}

bubble_aa = function(vec, radius_3d){
  
  return(which(vec <= radius_3d))
  
}


number_of_variants <- function(index,Uniprot_position){
  
  var_count.vec <- c()
  
  for(i in 1:length(index)){
    
    var_count.vec <- c(var_count.vec,variant_positions.vec %in% Uniprot_position[index[i] %>% unlist()] %>% 
                         which() %>% 
                         length()
    )
    
    
  }
  
  return(var_count.vec)
}


Patient_data.df <- read_delim("data/Patient_variants_SLC6A1_v4.txt", delim = "\t") %>% 
  select(-Transcript) %>% 
  mutate(AA_pos = as.numeric(AA_pos)) %>% 
  left_join(master.df %>% distinct(Transcript,Gene,AA_pos), by = c("AA_pos" = "AA_pos","Gene" = "Gene")) %>% 
  left_join(all_exchanges.df %>% distinct(Domain,Gene,AA_pos), by = c("AA_pos" = "AA_pos","Gene" = "Gene")) %>%  
  left_join(Functional_data.df %>% distinct(Gene,AA_pos,functional_effect), by = c("AA_pos" = "AA_pos","Gene" = "Gene"))  %>% 
  filter(selection_prio !=2) %>% 
  mutate(var = 1) %>%
  group_by(Aln_pos) %>% 
  summarise(path_count = sum(var))

Control_data.df <- read_delim("data/gnomad_variants.txt", delim = "\t") %>% 
  mutate(AA_ref = aaa(AA_ref),
         AA_alt = convert_aa(AA_alt)) %>% 
  left_join(master.df %>% distinct(Gene,AA_pos,Aln_pos), by = c("AA_pos" = "AA_pos","Gene" = "Gene")) %>% 
  mutate(var = 1) %>%
  group_by(Aln_pos) %>% 
  summarise(con_count = sum(var))

##enter the corresponding slc6a1 file here: SLC6A1_AF.txt in the pdb folder 
structure_coordiantes.df<- read_delim("data/pdb/6j8e_All_genes_structure_coordinates.txt", delim = "\t") %>% 
  filter(gene == "SLC6A1") %>% 
  left_join(master.df %>% filter(Gene == "SLC6A1") %>% distinct(Aln_pos,AA_pos), by = c("Uniprot_position" = "AA_pos")) %>% 
  select(Aln_pos,x,y,z) 
filter(!is.na(x))


Variant.df <- tibble(Aln_pos = 1:2181) %>% 
  left_join(Patient_data.df) %>% 
  left_join(Control_data.df) %>% 
  left_join(structure_coordiantes.df) %>% 
  filter(!is.na(x)) %>% 
  replace(is.na(.),0)

#select variants of interest 

#Input of the desired set of variants, here is just a dummy set of variants. We just need the Aminoacid Position of the variants (a Position can occur several times when more than one variant is located at the position) 
# variant_positions.vec <- read_csv("Grin_App/data_and_pdb/GRIN_db.csv") %>% 
#   filter(gene == "GRIN1",
#          !str_detect(variant_p,"dup")) %>% 
#   mutate(Uniprot_position = str_extract_all(variant_p,"[0-9]+") %>% unlist() %>%  as.numeric(),
#          variant = 1) %>% 
#   .$Uniprot_position  # Positions with selected variants 
dist.matrix <- Dist(Variant.df %>% select(x,y,z))

neighbour_patient <- c()
neighbour_gnomad <- c()
neighbour_p_outside <- c()
neighbour_g_outside <- c()

radius = 12

for(self.pos in 1:nrow(Variant.df)){ ### iterates each aminoacid position to calculate all distances 
  distance <- dist.matrix[self.pos,]
  
  neighbour_patient <- c(neighbour_patient,sum(Variant.df$path_count[which(distance<radius)]))
  neighbour_p_outside <- c(neighbour_p_outside,sum(Variant.df$path_count[which(distance>radius)]))
  
  neighbour_gnomad<- c(neighbour_gnomad,sum(Variant.df$con_count[which(distance<radius)]))
  neighbour_g_outside<- c(neighbour_g_outside,sum(Variant.df$con_count[which(distance>radius)]))
  
}


output.df <- tibble(p_in = neighbour_patient , p_out = neighbour_p_outside , g_in = neighbour_gnomad, g_out = neighbour_g_outside, Aln_pos = Variant.df$Aln_pos)


total_number_patient <- sum(output.df$p_in,output.df$p_out)
total_number_gnomad <- sum(output.df$g_in,output.df$g_out)
fisher_odds <- c()
fisher_pvalue <- c()


for(pos in 1:nrow(output.df)){
  fisher_input.matrix <- matrix(c(output.df$p_in[pos],output.df$g_in[pos],
                                  output.df$p_out[pos],
                                  output.df$g_out[pos]), byrow = T, ncol = 2) ##creates fisher matrix which looks like ordered by row: 
  #Number of patients in the sphere; Number of Gnomad variants in the sphere
  #Number of patients outside the sphere; Number of Gnomad variants outside the sphere
  
  fishtest <- fisher.test(fisher_input.matrix)  ##perfrom fisher test 
  fisher_odds <- c(fisher_odds,as.numeric(fishtest$estimate))  ##save odds ratios in a vector 
  fisher_pvalue <- c(fisher_pvalue,as.numeric(fishtest$p.value))
  
}

output.df %>% 
  mutate(pvalue = fisher_pvalue,
         odds = fisher_odds,
         PER3D = ifelse(odds >1 & pvalue < 0.05,"PER3D","No-PER3D")) %>% 
  write_delim("data/per3d.txt",delim = "\t")