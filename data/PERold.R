#For one Gene 

Patient_data.df <- read_delim("data/Patient_variants_SLC6A1_v5.txt", delim = "\t")

Control_data.df <- read_delim("data/gnomad_variants.txt", delim = "\t")

patient.df <- Patient_data.df %>% 
  group_by(AA_pos,Gene,AA_alt) %>% 
  # filter(Gene %in% c("SCN1A","SCN2A","SCN3A","SCN8A")) %>% 
  summarise(pscount = n()) %>% 
  group_by(AA_pos,Gene) %>% 
  summarise(pscount = sum(pscount))

control.df <- Control_data.df %>% 
  group_by(AA_pos,Gene,AA_alt) %>% 
  # filter(Gene %in% c("SCN1A","SCN2A","SCN3A","SCN8A")) %>% 
  summarise(gscount = n()) %>% 
  group_by(AA_pos,Gene) %>% 
  summarise(gscount = sum(gscount))


aln_master.df <- read_delim("data/master_table.txt", delim = "\t") %>% 
  filter(Gene %in% c("SLC6A1")) %>% 
  distinct(AA_pos,Gene) %>% 
  arrange(AA_pos,Gene) %>% 
  mutate(spec_aln_pos = 1:nrow(.)) %>% ##alignment position of only scn1a,2,3,8 , ignoring the other ones
  left_join(patient.df) %>% 
  left_join(control.df) %>% 
  select(AA_pos,spec_aln_pos,pscount,gscount,Gene) %>% 
  replace(is.na(.),0) %>% 
  group_by(Gene) %>% 
  mutate(in_ps = sliding_window_2d(4,pscount),
         in_gs = sliding_window_2d(4,gscount)) %>% 
  ungroup() %>% 
  group_by(Gene) %>% 
  mutate(out_ps = sum(pscount) -in_ps,
         out_gs = sum(gscount) -in_gs,
         pvalue = 999,
         odds = 999,
         selected_pos = 0) %>% 
  ungroup()

sliding_window_2d <- function(windowsize,vec){
  #windowsize = number of aminoacids to be considered at each side - windowsize of 3 covers 7 aminoacids in total !
  
  window_sum <- c()
  
  for(i in 1:length(vec)){
    
    
    window_sum <- c(window_sum,
                    sum(vec[ifelse(i-windowsize<1,1,i-windowsize):ifelse(i+windowsize>length(vec),length(vec),i+windowsize)]))
    
  }
  
  return(window_sum) # the output will be the number of variants that are inside the window  for each position
  
}

#aln_master.df <- tibble()

for(i in 1:nrow(aln_master.df)){
  
  if(i%%5 == 0 | i == 1){
    mat <- matrix(c(aln_master.df$in_ps[i],aln_master.df$in_gs[i],aln_master.df$out_ps[i],aln_master.df$out_gs[i]), ncol = 2)
    
    pvalue <- fisher.test(mat) %>% .$p.value
    odds<-fisher.test(mat) %>% .$estimate
    
    aln_master.df$selected_pos[i] <- 1
  }
  
  
  aln_master.df$pvalue[i] <- pvalue
  aln_master.df$odds[i] <- odds
}

# for(gene_sel in aln_master.df$Gene %>% unique()){
#   
#   aln_master_gene.df <- aln_master.df %>% filter(Gene == gene_sel)
#   
#   for(i in 1:nrow(aln_master_gene.df)){
#     
#     mat <- matrix(c(aln_master_gene.df$in_ps[i],aln_master_gene.df$in_gs[i],aln_master_gene.df$out_ps[i],aln_master_gene.df$out_gs[i]), ncol = 2)
#     
#     pvalue <- fisher.test(mat) %>% .$p.value
#     odds<-fisher.test(mat) %>% .$estimate
#     
#     aln_master_gene.df$pvalue[i] <- pvalue
#     aln_master_gene.df$odds[i] <- odds
#     
#     aln_master.df$selected_pos[i] <- 0
#     
#     if(i%%5 == 0 | i == 1){
#       
#       print(i)
#       
#       aln_master.df$selected_pos[i] <- 1
#     }
#   }
#   
#   
#   final_master.df <- rbind(final_master.df,aln_master_gene.df)
#   
#   
# }


#adapted for SLC6A1 or any portal with just one gene. Ensure that final_master.df$Gene %>% unique() shows only one gene 

per_family2d.df <- aln_master.df %>% 
  mutate(pvalue= ifelse(pvalue* sum(selected_pos) > 1,1,pvalue* sum(selected_pos)),
    odds = ifelse(odds == Inf,max(odds[odds != Inf]),odds),
    per = ifelse(pvalue <0.05 & odds >1,"PER","neutral")) 


for(i in 1:nrow(per_family2d.df)){
  
  if(per_family2d.df$selected_pos[i] == 1 & per_family2d.df$pvalue[i] <0.05 & per_family2d.df$odds[i]){
    
    index <- c(ifelse((i-4)<1,1,i-4):ifelse((i+4)>nrow(per_family2d.df),nrow(per_family2d.df),i+4))
    
    per_family2d.df$per[index] <- "PER"
    
  }
  
}


write_delim(per_family2d.df, "data/per2d_genewise.txt", delim = "\t")
