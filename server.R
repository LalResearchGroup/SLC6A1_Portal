################################################
################# SLC6A1 PORTAL ##################
################################################
################## LIBRARIES ################

library(shiny)
library(plotly)
library(readxl)
library(readr)
library(r3dmol)
library(shinyWidgets)
library(shinydashboard)
library(shinydashboardPlus)
library(DT)
library(tidyverse)
library(seqinr)
library(bio3d)
library(RColorBrewer)
library(plyr)
library(tippy)
library(vembedr)
library(shinyhelper)
library(rsconnect)
library(pacman)
library(Rfast)

############## FUNCTIONS, STYLE ############

numextract <- function(string) {
    as.numeric(str_extract(string, "\\-*\\d+\\.*\\d*"))
}

convert_aa <- function(aa){
  
  return(ifelse(aa == "*","Stop",aaa(aa)))
  
}

log_fun <- function(x) { 
    x = case_when(
        as.numeric(x) > 0 ~ log10(x), #log10 transformation to make the visuals clean 
        as.numeric(x) == 0 ~ 0,
        as.numeric(x) < 0 ~ log10(x * -1) * -1, 
        TRUE ~ -9999)
    return(x)
    
}

plotly_font <- list(
    family = "sans-serif",
    size = 15)

goodbye <- c("zoomIn2d", "zoomOut2d", "hoverClosestGeo", "hoverClosestGl2d",
             "hoverClosestCartesian","lasso2d","select2d","resetScale2d",
             "hoverCompareCartesian", "hoverClosestPie","toggleSpikelines")

# coloring for colorblindness 
lolliplot_fill_scheme <-  #("Missense"="#D55E00","PTV"="#0072B2","Control" ="#000000", "PER" = "red"))+
    c(
        "N-Terminal Cytoplasmic" = "#FFFF00",
        "Transmembrane" = "#66FFFF",
        "Linker" = "lightgray",
        "C-Terminal-Cytoplasmic"="#FFB266",
        "no" = "#FFFFFF",
        "Missense" = "#D55E00",
        "Synonymous" = "#D55E00",
        "PTV" = "#0072B2",
        "Other" = "#0072B2"
        )



######Functions######

add_odds <- function(g1_in,g1_out,g2_in,g2_out,label_sel){
  
  out.df <- tibble()
  
  for(i in 1:length(g1_in)){
    
    fisher_out <- matrix(c(g1_in[i],g1_out[i],g2_in[i],g2_out[i]), ncol = 2) %>% fisher.test()
    
    out.df <-  rbind(out.df,
                     tibble(odds = fisher_out$estimate, pvalue = fisher_out$p.value, lCI = fisher_out$conf.int[1], uCI = fisher_out$conf.int[2], label = label_sel[i]))
    
    
  }
  
  # out.df <- out.df %>%
  #   mutate(pvalue_adj = pvalue*nrow(.))
  
  return(out.df)
  
}


#Basic Information 

Phenotype_fac_1.fun <- function(select_gene,phenotype,color_sel){
  
  plot_ly(Patient_data.df %>% 
            filter(Gene == select_gene) %>% 
            dplyr::rename(phenotype_sel = phenotype) %>% 
            select(phenotype_sel) %>%
            arrange(phenotype_sel) %>% 
            mutate(phenotype_sel = factor(phenotype_sel, levels = unique(phenotype_sel))) %>% 
            filter(!is.na(phenotype)) %>% 
            dplyr::group_by(phenotype_sel) %>% 
            dplyr::summarise(n = n()) %>% 
            assign("save",.,envir = .GlobalEnv), 
          values=~n, labels=~phenotype_sel, type="pie",
          textposition="inside", 
          #colors = "Set2",
          textinfo="label+percent",
          insidetextfont = list(color = '#333333'),
          hoverinfo = "text",
          text= ~ paste(n, "individuals"),
          marker=list(colors= basic_phenotype_colors[1:nrow(save)],
                      line = list(color = '#FFFFFF', width = 1)), showlegend = FALSE) %>% 
    layout(title=paste(""),font=plotly_font) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  
}


Phenotype_fac_2.fun <- function(select_gene,phenotype,color_sel){
  
  plot_ly(Patient_data.df %>% 
            dplyr::rename(phenotype_fac = phenotype) %>% 
            filter(phenotype_fac != "NA") %>% 
            mutate(phenotype_fac =ifelse(phenotype_fac == "Yes"," Yes",phenotype_fac)) %>% 
            select(phenotype_fac) %>%
            arrange(phenotype_fac) %>% 
            mutate(phenotype_sel = factor(phenotype_fac, 
                                          levels = unique(phenotype_fac))) %>% # c("Severe DD/ID", "Mild DD/ID", "learning disability", "Unclassified DD", "Normal", "Moderate DD/ID"),
                                          # ordered = TRUE)) %>% 
            dplyr::group_by(phenotype_fac) %>%
            dplyr::summarise(n = n()) %>% 
            mutate(n_gene = sum(n))%>% 
            assign("save",.,envir = .GlobalEnv),
          x = ~ phenotype_fac, 
          y = ~ round(n, digits = 2), 
          color = ~ phenotype_fac, 
          colors = color_sel,
          type = "bar", 
          hoverinfo = "text", showlegend = FALSE,
          text= ~ paste0(round(n, digits = 2), " (" ,n," individuals)")) %>% 
    layout(title="", 
           font=plotly_font,
           xaxis = list(title="",showline = T, tickangle = 45),
           yaxis = list(title="N of individuals",showline = T),
           margin = list(b = 160)) %>%
    config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  
  
}

Phenotype_fac_3.fun <- function(select_gene,phenotype,color_sel,filter_var){
  
  plot_ly(Patient_data.df %>% 
            dplyr::rename(phenotype_fac = phenotype) %>% 
            filter(Vartype == filter_var) %>% 
            select(phenotype_fac) %>%
            dplyr::arrange(phenotype_fac) %>% 
            dplyr::mutate(phenotype_sel = factor(phenotype_fac, levels = unique(phenotype_fac))) %>% 
            dplyr::group_by(phenotype_fac) %>% 
            dplyr::summarise(n = n()) %>% 
            dplyr::mutate(n_gene = sum(n))%>% 
            assign("save",.,envir = .GlobalEnv),
          x = ~ phenotype_fac, 
          y = ~ round(n, digits = 2), 
          color = color_sel,
          colors = color_sel,
          type = "bar", 
          hoverinfo = "text", showlegend = FALSE,
          text= ~ paste0(round(n, digits = 2), " (" ,n," individuals)")) %>% 
  layout(title = paste0(filter_var," variants"), 
         font=plotly_font,
         xaxis = list(title="",showline = T, tickangle = 45),
         yaxis = list(title="N of individuals",showline = T),
         margin = list(b = 160, t = 50)) %>%
    config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE) 
  
}

## Seizure onset

Onset_days_1.fun <- function(select_gene,phenotype,colors_sel){
  
  plot_ly(Patient_data.df %>% 
            dplyr::rename(Onset_days = phenotype) %>% 
            dplyr::mutate(p_variant = paste0("p.",AA_ref,AA_pos,AA_alt)) %>% 
            filter(Vartype != "Synonymous",
                   Vartype != "INDEL",
                   Vartype != "Splice",
                   Vartype != "CNV",
                   Vartype != "Complex") %>%
            dplyr::mutate(ID = paste0(ifelse(Vartype == "Missense","Missense","PTV"))) %>% 
            filter(!is.na(Onset_days)) %>% 
            dplyr::group_by(ID) %>% 
            dplyr::mutate(counts = n()) %>% 
            dplyr::ungroup() %>% 
            dplyr::mutate(ID_save = ID,
                   ID = paste0(ID,"\n(N = ",counts, ")")) %>% 
            assign("save",.,envir = .GlobalEnv) %>% 
            dplyr::arrange(ID_save) %>% 
            dplyr::mutate(ID = factor(ID, levels = unique(ID))),
          y = ~Onset_days, color = ~ID, colors = colors_sel, type = "box",
          boxpoints = "all", jitter = 0.3,
          pointpos = 0, hoverinfo = "text", showlegend = FALSE,
          text= ~paste0(round(Onset_days, digits = 2), " Days, ", Protein)) %>% 
    layout(font=plotly_font,  
           title="",
           xaxis = list(title="", tickangle = 45, showline = T),
           yaxis = list(type = "log",
                        title = "Seizure onset (months)",
                        tickvals = list(0.01,0.5,1,5,10,20,50,100,200,500,1000,2000,5000,10000,20000,50000),
                        tickmode = "array",
                        showline = T
           )) %>%
    config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  
}


Onset_days_2.fun <- function(select_gene,phenotype,colors_sel){
  plot_ly(Patient_data.df %>% 
            filter(Gene == select_gene,
                   !is.na(Sz_onset)) %>% 
            dplyr::rename(phenotype_fac = phenotype) %>% 
            filter(Vartype != "Synonymous",
                   Vartype != "Complex",
                   Vartype != "Missense mosaic") %>% 
            dplyr::mutate(ID = paste0(phenotype_fac,ifelse(Vartype == "Missense","Missense","PTV"))) %>% 
            dplyr::arrange(ID) %>% 
            dplyr::mutate(ID = factor(ID, levels = unique(ID))) %>% 
            dplyr::group_by(ID) %>% 
            dplyr::summarise(n = n()) %>% 
            dplyr::mutate(n_gene = sum(n))%>% 
            assign("save",.,envir = .GlobalEnv),
          x = ~ ID, 
          y = ~ n, 
          color = ~ ID, 
          colors = colors_sel,
          type = "bar", 
          hoverinfo = "text", showlegend = FALSE,
          text= ~ paste0(round(n/n_gene, digits = 2), " (" ,n," individuals)")) %>% 
    layout(title="", 
           font=plotly_font,
           xaxis = list(title="",tickangle=45,showline = T),
           yaxis = list(title="Number of patients",showline = T)) %>%
    config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  
}

basic_onset_legend <- function(){
  legend <- data.frame(x=c(1,3), y=c(3,3), text=c("Missense", "PTV       "))
  plot <- ggplot(legend, aes(x=x, y=y, color=text))+
    geom_point(size = 8)+
    scale_color_manual(values = c("Missense"="#BEBADA","PTV       "="#FDB462"))+
    ylim(c(2.9,3.1))+
    xlim(c(0.8,5.2))+
    theme_void()+
    geom_text(aes(label=text), hjust=-0.2, color="black", size =5)+
    theme(legend.position = "none")
  
  
  return(plot)
}

#variant analysis< 
extract_gnomad_features <- function(Control_data.df,selected.df,variable,label){

  if(label == "exchange"){
    
    control_int.df <- Control_data.df %>% filter(Gene == selected.df$Gene[1], AA_pos == selected.df$AA_pos[1], AA_alt == selected.df$AA_alt[1])
    
    control_int2.df <- Control_data.df %>% filter(Gene == selected.df$Gene[1], AA_pos == selected.df$AA_pos[1]) 
    
    if(variable == "Allele count"){
      
      out <- ifelse(nrow(control_int.df) == 0,0,control_int.df$Allele_count)
      
    }else if(variable == "Allele freq"){
      
      out <- ifelse(nrow(control_int.df) == 0,0,control_int.df$Allele_freq)
      
    }
  }else if(label == "position"){
    
    control_int.df <- Control_data.df %>% filter(Gene == selected.df$Gene[1], AA_pos == selected.df$AA_pos[1]) 
    
    if(variable == "Allele count"){
      
      out <- ifelse(nrow(control_int.df) == 0,0,sum(control_int.df$Allele_count))
      
    }else if(variable == "Allele freq"){
      
      out <- ifelse(nrow(control_int.df) == 0,0,sum(control_int.df$Allele_freq))
      
    }
    
  }
  
  return(out)
  
  
}

##map paraz score 
map_paraz <- function(data_scores.df){
  plot <- ggplot(data_scores.df , aes(x = AA_pos, y = Paraz_score, fill = group)) +
    geom_bar(stat = "identity") +
    theme_bw() + theme(panel.border = element_blank(),
                       legend.position = "none",
                       legend.title = element_blank()) +
    scale_fill_manual(values = c("gray","indianred"))+
    coord_cartesian(ylim = c((min(data_scores.df$Paraz_score,na.rm = T)-0.25),
                             max(data_scores.df$Paraz_score,na.rm = T)+0.25),
                    expand = FALSE) +
    labs(title="Paralog conservation",
         y="Parazscore",
         x = paste0("Amino acid sequence"))
  
  return(plot)
}


##map mtrscore 
map_mtr <- function(data_scores.df){
  
  mtr_threshold = data_scores.df$MTR_score %>% sort() %>% .[round(length(.)/4)] #lowest 25% of values. This threshold is dervied from the publication of the score 
  mtr_threshold_hard = data_scores.df$MTR_score %>% sort() %>% .[round(length(.)/20)] #lowest 5% of values. This threshold is dervied from the publication of the score 
  
  plot <- ggplot(data_scores.df, aes(x = AA_pos, y = MTR_score, colour=(MTR_score<mtr_threshold))) +
    #geom_line(aes(group=group)) +
    geom_point(aes(x = AA_pos2), color = "purple")+
    #scale_color_manual(values = c("gray","indianred1"),
    #                   labels = c("Tolerant region", "Intolerant region")) +
    geom_hline(aes(yintercept=mtr_threshold), colour="orange", linetype="dashed")  +
    geom_hline(aes(yintercept=mtr_threshold_hard), colour="indianred1", linetype="dashed")  +
    coord_cartesian(expand = FALSE) +
    theme_bw() + theme(panel.border = element_blank(),
                       legend.position = "right",
                       legend.direction = "vertical") +
    labs(title="Missense Tolerance Ratio",
         y="MTR",
         x = paste0("Amino acid sequence"),
         colour = "Patient Variants")
  
  return(plot)
}



#Research 3d mapping
map_var_3d <- function(data,Gene_sel,gnomad_bool,pdb_sel,structure_coordinates,sub_color_i){

  variant.df <- data %>%
    filter(Gene == Gene_sel) %>%
    filter(Vartype == "Missense") %>%
    mutate(label = "pathogenic") %>%
    dplyr::group_by(AA_pos,AA_ref,Gene) %>%
    dplyr::summarise(n_occ = n()) %>%
    select(AA_pos,AA_ref,n_occ,Gene)
  
  gnomad.df <- Control_data.df %>%
    filter(Gene == Gene_sel) %>%
    dplyr::group_by(AA_pos,AA_ref,Gene) %>%
    dplyr::summarise(n_occ = n()) %>%
    select(AA_pos,AA_ref,n_occ,Gene)
  
  structure.df <- read_delim(structure_coordinates,delim = "\t") %>%
    mutate(Aminoacid = aaa(Aminoacid)) %>%
    select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
  
  
  variant.df <- variant.df %>%
    left_join(structure.df,by = c("AA_pos" = "Uniprot_position","AA_ref" = "Aminoacid","Gene" = "gene")) %>%
    mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not
    filter(struc_cov == "yes")
  
  gnomad.df <- gnomad.df %>%
    left_join(structure.df,by = c("AA_pos" = "Uniprot_position","AA_ref" = "Aminoacid","Gene" = "gene")) %>%
    mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not
    filter(struc_cov == "yes")
  
  sub_color <- c("#e2f970","#6fbbf7","#ee6c71","#ffbc5a","#bf73cc")
  sub_scale <- c(1.2,0.8)
  struc_color <- "wheat"
  
  rot = 270
  rot_axis = "x"
  spin_axis = "vy"
  
  #Specify yourself- color of the cartoon per subunit
  subunit_color <- c("wheat","white") #
  
  print(variant.df)
  print("now gnomad")
  print(gnomad.df)
  #Model for the protein complex
  
  modelo <- r3dmol(
    viewer_spec = m_viewer_spec(
      cartoonQuality = 10,
      lowerZoomLimit = 50,
      upperZoomLimit = 1000
    )
  )
  
  modelo <- modelo %>% m_add_model(data = pdb_sel, format = "pdb")
  
  # Zoom to encompass the whole scene
  modelo <- modelo %>% m_zoom_to() %>%
    # Set color o cartoon representation
    m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure
    # Set subunit colors
    m_set_style(
      sel = m_sel(chain = c("A")),
      style = m_style_cartoon(color = subunit_color[2])
    ) %>%
    # visualize variants grin1
    m_set_style(
      sel = m_sel(resi = variant.df$Position_in_structure,
                  atom = "CA",
                  chain = c("A")),
      style = m_style_sphere(colorScheme = NULL,
                             color = sub_color[sub_color_i],
                             scale = sub_scale[1])
    )
  
  if  (gnomad_bool == TRUE) {
    
    modelo <- modelo %>% m_set_style(
      sel = m_sel(resi = gnomad.df$Position_in_structure,
                  atom = "CA",
                  chain = c("A")),
      style = m_style_sphere(colorScheme = NULL,
                             color = "#333333",
                             scale = sub_scale[2]))
    
  }

  return(modelo) 
}     

#functional data #####
func_dist <- function(data.df,chosen_dist,title_sel,y_sel,color_sel){
  
  plot <- plot_ly(data.df %>% 
                    filter(Vartype == "Missense") %>% 
                    mutate(act_group = case_when(`GABA uptake (vs wt)`<10~"<10",
                                                 `GABA uptake (vs wt)`<42.8~"10-42.8",
                                                 `GABA uptake (vs wt)`>42.8~">42.8")) %>% 
                    mutate(act_group = factor(act_group, levels = c("<10","10-42.8",">42.8"))) %>% 
                    #act_group = factor(act_group, levels = c("<10","10-42.8",">42.8"))) %>% 
                    dplyr::rename(dist = chosen_dist) ,
                  x = ~act_group, y = ~dist,
                  type = "box",
                  boxpoints = "all",
                  jitter = 0.2,
                  pointpos = 0, 
                  fillcolor = color_sel,
                  marker = list(color = "black"),
                  line = list(color = "black")) %>% 
    #fillcolor = "khaki") %>% 
    layout(title=title_sel, 
           font=plotly_font,
           xaxis = list(title="Average GABA uptake rate in relation to WT ( %)"),
           margin = list(t = 50),
           yaxis = list(title=y_sel)) %>% 
  config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)

return(plot)

}


map_func_var <- function(data.df,group_sel,color_sel){
  
  structure.df <- read_delim("data/pdb/7sk2_struc.txt",delim = "\t") %>%
    mutate(Aminoacid = aaa(Aminoacid)) %>%
    select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
  
  
  var_map.df <- data.df %>% 
    mutate(act_group = case_when(`GABA uptake (vs wt)`<10~"<10",
                                 `GABA uptake (vs wt)`<42.8~"10-42.8",
                                 `GABA uptake (vs wt)`>42.8~">42.8")) %>% 
    distinct(act_group,AA_pos,AA_ref) %>% 
    left_join(structure.df,by = c("AA_pos" = "Uniprot_position","AA_ref" = "Aminoacid"))%>%
    filter(act_group== group_sel) %>% 
    filter(!is.na(Position_in_structure)) %>% 
    distinct(Position_in_structure)
  
  
  sub_color <- c("<10" = "red",
                 "10-42.8" = "orange",
                 ">42.8" = "yellow")
  
  sub_scale <- c(1,1)
  struc_color <- "white"
  
  rot = 270
  rot_axis = "x"
  spin_axis = "vy"
  
  #Specify yourself- color of the cartoon per subunit
  subunit_color <- c("wheat","white") #first color for GRIN1 second or GRIN2A
  # 
  # #Model for the protein complex
  
  modelo <- r3dmol(
    viewer_spec = m_viewer_spec(
      cartoonQuality = 10,
      lowerZoomLimit = 50,
      upperZoomLimit = 1000
    )
  )
  
  modelo <- modelo %>% m_add_model(data = "data/pdb/7sk2.pdb", format = "pdb")
  
  # Zoom to encompass the whole scene
  modelo <-modelo %>% m_zoom_to() %>%
    # Set color o cartoon representation
    m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure
    # Set subunit colors
    m_set_style(
      sel = m_sel(chain = c("A")),
      style = m_style_cartoon(color = subunit_color[2])
    ) %>%
    # visualize variants all
    m_set_style(
      sel = m_sel(resi = var_map.df %>% .$Position_in_structure,
                  atom = "CA",
                  chain = c("A")),
      style = m_style_sphere(colorScheme = NULL,
                             color = color_sel,
                             scale = sub_scale[1])
    )
  
  return(modelo)
}

## substitute for seqinr function
three_to_one_aa <- function(sequence){
  
  sapply(sequence, function(x){
    
    ifelse(is.na(x),NA,seqinr::aaa(x))
    
  }) %>% 
    as.vector() %>% 
    return()
  
}


############## DATA ############



#Load domain data 
Domain_data.df <- read_delim("data/Domain_cornelius.txt", delim = "\t") %>% 
  dplyr::rename(Domain = Cornelius_level1, Domain_color = Cornelius_level2) %>% 
  dplyr::select(Aln_pos,Domain,Domain_color)

#Load all possible exchanges 
all_exchanges.df <- read_delim("data/master_table_exchanges.txt",delim = "\t") %>% 
  left_join(Domain_data.df %>% distinct(Aln_pos,Domain,Domain_color)) %>% 
  mutate(AA_ref = aaa(AA_ref),
         AA_alt = convert_aa(AA_alt))

#Load master table 
master.df <- read_delim("data/master_table.txt", delim = "\t") 

#load domains
Domain_gene.df <- master.df %>% 
  left_join(Domain_data.df %>% distinct(Aln_pos,Domain,Domain_color)) %>% 
  distinct(Gene,AA_pos,Domain,Domain_color) %>% 
  ungroup() %>% 
  mutate(Domain = Domain_color)


#PER-2D
per_family2d.df <- read_delim("data/per2d_genewise.txt", delim = "\t")

per2d_for_var_analysis.df <- read_delim("data/per2d_genewise.txt", delim = "\t")

per2d.df <- read_delim("data/per2d_genewise.txt", delim = "\t") %>% 
  mutate(Hotzone_2D = ifelse(per == "PER","PER","No-PER")) %>% 
  select(Hotzone_2D,AA_pos,pvalue,odds) %>% 
  dplyr::rename(pvalue_per2d = "pvalue",
         odds_per2d = "odds")
#PER-3D
per3d.df <- read_delim("data/per3d.txt", delim = "\t") %>% 
  select(PER3D,AA_pos,pvalue,odds) %>% 
  dplyr::rename(pvalue_per3d  = "pvalue",
         odds_per3d  = "odds")

# Adjust all_exchanges.df
exchanges.df <- all_exchanges.df %>% 
  mutate(AA_ref = a(AA_ref),
         AA_alt = a(AA_alt)) %>% 
  select(AA_pos, AA_ref, AA_alt, Vartype)


dist_imp_features.df <- read_delim("data/Distance_imp_features.txt",delim = "\t") %>% 
  mutate(Vartype = "Missense") %>% 
  mutate(AA_ref = aaa(AA_ref))

#Load functional data mermer paper
biomarin.df <- read_csv("data/Functional_data_biomarin.csv") %>% 
  dplyr::rename(cDNA_ref ="Ref",
                cDNA_alt = "Alt",
                cDNA_pos = "Pos",
                Vartype = "variant_impact",
                `GABA uptake (vs wt)` = "Avg_PercentWT") %>% 
  select(Gene,cDNA_pos,cDNA_ref,cDNA_alt,AA_pos,AA_ref,AA_alt,Vartype,`GABA uptake (vs wt)`) %>%
  mutate(Vartype = case_when(Vartype == "missense" ~"Missense",
                             Vartype == "frameshift" ~"PTV",
                             Vartype == "stop gained" ~"PTV",
                             Vartype == "inframe indel" ~"Indel",
                             Vartype == "synonymous" ~"Synonymous")) %>% 
  mutate(
    var_id = paste0(AA_ref,AA_pos,AA_alt))


mermer.df <- read_delim("data/Functional_data_mermer.txt", delim = "\t") %>%
  mutate(Vartype = ifelse(AA_alt != "X","Missense","PTV"),
         cDNA_alt = NA,
         AA_ref = aaa(AA_ref),
         AA_alt = aaa(AA_alt)) %>% 
  select(Gene,cDNA_pos,cDNA_ref,cDNA_alt,AA_pos,AA_ref,AA_alt,Vartype,`GABA uptake (vs wt)`,
         `Surface expression (vs wt)`,
         `Total expression (vs wt)`,
         `Relative uptake/surface expression`,
         `Relative surface expression/total expression`) %>% 
  mutate(var_id = paste0(AA_ref,AA_pos,AA_alt)) %>% 
  mutate(`GABA uptake (vs wt)` = `GABA uptake (vs wt)`*100)


functional.df <- biomarin.df %>% 
  rbind(mermer.df %>% filter(!(var_id %in% biomarin.df$var_id)) %>% 
          select(-`Surface expression (vs wt)`, -`Total expression (vs wt)`,-`Relative uptake/surface expression`,-`Relative surface expression/total expression`)) %>% 
  left_join(mermer.df %>% select(var_id, `Surface expression (vs wt)`, `Total expression (vs wt)`,`Relative uptake/surface expression`,`Relative surface expression/total expression`), by = c("var_id" = "var_id")) %>% 
  select(AA_pos,AA_ref,AA_alt,Vartype,`Surface expression (vs wt)`, `Total expression (vs wt)`,`Relative uptake/surface expression`,`Relative surface expression/total expression`,`GABA uptake (vs wt)`)


functional_only.df <- functional.df%>% 
  left_join(Domain_gene.df %>% distinct(Domain,Gene,AA_pos,Domain_color), by = c("AA_pos" = "AA_pos")) %>%
  left_join(per3d.df) %>% 
  left_join(per2d.df) %>% 
  left_join(dist_imp_features.df)

# Exchanges for patient data
exchanges <- all_exchanges.df %>% 
  select(cDNA_pos, cDNA_ref, cDNA_alt, AA_pos, AA_ref, AA_alt, Vartype)

#Load patient and control data 
Patient_data.df <- read_delim("data/Patient_variants_SLC6A1_v8.txt", delim = "\t") %>% 
  select(-Transcript) %>% 
  mutate(AA_pos = as.numeric(AA_pos)) %>% 
  ##specific to each dataset
  dplyr::rename(Sz_onset = "Age at seizure onset (months)",
                Epilepsy = "Epilepsy") %>% 
  left_join(master.df %>% distinct(Transcript,Gene,AA_pos), by = c("AA_pos" = "AA_pos","Gene" = "Gene")) %>% 
  left_join(Domain_gene.df %>% distinct(Domain,Gene,AA_pos,Domain_color), by = c("AA_pos" = "AA_pos","Gene" = "Gene")) %>%
  left_join(per3d.df) %>% 
  left_join(per2d.df) %>% 
  mutate(AA_ref = ifelse(!is.na(AA_ref),AA_ref,"XXX") %>% aaa(), ##warnings due to none matching aminoacids are fine 
         AA_alt = ifelse(!is.na(AA_alt),AA_alt,"XXX") %>% aaa(),
         #AA_alt_complex = ifelse(Vartype == "Missense",AA_alt,AA_alt_complex),
         cDNA = ifelse(!is.na(cDNA_pos), paste0("c.",cDNA_pos,cDNA_ref,">",cDNA_alt), "Not available"),
         #Protein = paste0("p.",AA_ref,AA_pos,AA_alt_complex),
         Protein = Original_AA_change,
         Inheritance = ifelse(is.na(Inheritance),"NA",Inheritance),
         Epilepsy = ifelse(Epilepsy %in% c("No","Yes"),Epilepsy,"NA")) %>% 
  dplyr::rename(Autism = "Autistic traits",
                Epilepsy_syndrome = "Epilepsy Syndrome Classification",
                ID_after_sz_onset = "Cognitive Level AFTER Seizure Onset") %>% 
  mutate(Autism = ifelse(is.na(Autism),NA,Autism)) %>% 
  left_join(functional.df %>% filter(Vartype == "Missense"), by = c("AA_pos" = "AA_pos","AA_alt" = "AA_alt","AA_ref" = "AA_ref","Vartype" = "Vartype")) %>% 
  left_join(dist_imp_features.df) %>% 
  left_join(exchanges, by = c("cDNA_pos", "cDNA_ref", "cDNA_alt")) %>% 
  mutate(Vartype = case_when(Vartype.x == "CNV" ~ "CNV",
                             Vartype.x == "INDEL" ~ "INDEL",
                             Vartype.x == "PTV" ~ "PTV",
                             Vartype.x == "?" ~ "Missense",
                             Vartype.x == "Missense" ~ "Missense",
                             Vartype.x == "Splice" ~ "Splice")) %>% 
  select(-Vartype.x, -AA_pos.y, -AA_ref.y, -AA_alt.y, -Vartype.y) %>% 
  dplyr::rename(AA_pos = AA_pos.x, AA_ref = AA_ref.x, AA_alt = AA_alt.x)

Patient_data_missense_only.df <- Patient_data.df %>% 
  filter(Vartype == "Missense")

# All functionally tested variants with clinical data if available
functional.df.missense_patient.df <- functional.df.missense %>% 
  left_join(Patient_data.df) %>% 
  distinct(AA_pos, AA_ref, AA_alt, `GABA uptake (vs wt)`, .keep_all = T)

Control_data.df <- read_delim("data/gnomad_variants.txt", delim = "\t") %>% 
  mutate(AA_ref = aaa(AA_ref),
         AA_alt = convert_aa(AA_alt)) %>%
  left_join(per3d.df) %>% 
  left_join(per2d.df) %>% 
  left_join(Domain_gene.df %>% distinct(Domain,Gene,AA_pos,Domain_color), by = c("AA_pos" = "AA_pos","Gene" = "Gene")) 

#Load Scores
#Paraz/MTR
paraz_mtr.df <- read_delim("data/mtr_paraz_slc6a1.txt",delim = "\t")

#hotzones3D on structure 
#PER3D_struc.df <- read_delim("data/pdb/6j8e_varburden.txt", delim = "\t")

#Load ClinVar data for Variant Analysis
clinvar.df <- read_delim("data/Clinvar_links_SLC6A1.txt", delim = "\t")

##### Variables #####
  #Basic Information
  basic_gene1 = "SLC6A1"
  basic_phenotype_fac = "Epilepsy"
  basic_phenotype_num = "Sz_onset"
  
  
  phenotype_name1 <- "Epilepsy"
  phenotype_name2 <- "Seizure onset (months)"
  phenotype_name3 <- "Seizure onset (months)"
  
  basic_phenotype_colors <- RColorBrewer::brewer.pal(20,"Set3") # Warning occurs at the moment, can be ignored for now 
  
  basic_phenotype_colors_ID <- c("black","grey","#994714")
  
  basic_phenotype_colors_epi <- c("No" = "grey",
                              " Yes" = "purple",
                              "Yes" = "purple")
  basic_phenotype_colors_autism <- c("No" = "grey",
                                  " Yes" = "#EFC56F",
                                  "Yes" = "#EFC56F")
  colors_gene1 <-  c("#BEBADA","#FDB462","#BEBADA","#FDB462")
 
##### Variant Analysis variable #####
  variant_title1 <- "Epilepsy"
  variant_title2 <- "Autism"
  variant_title3 <- "Seizure onset"

#####Research variable #####
  Gene_colors <-  c("Patient"="red", 
                    "Control" = "blue",
                    "Other" = "#333333")
  
  research_phenotype1_title <- "Number of patient variants per unit"
  research_phenotype2_title <- "Epilepsy"
  research_phenotype3_title <- "Seizure onset"
  research_phenotype4_title <- "Autism"
  research_phenotype5_title <- "Cognitive development"
  research_phenotype6_title <- "Epilepsy syndrome classification"
  
  research_functional1_title <- "Variants with molecular function assessment"
  #3d_mapping Genotype interface 
  pdb_sel_gene1 = "data/pdb/SCN1A_model.pdb1"
  structure_coordinates_gene1 <- "data/pdb/7dtd_structure_coordinates.txt"
  
  pdb_sel_gene2 = "data/pdb/SCN2A_model.pdb1"
  structure_coordinates_gene2 <- "data/pdb/6j8e_structure_coordinates.txt"
  
  pdb_sel_gene3 = "data/pdb/SCN2A_model.pdb1"
  structure_coordinates_gene3 <- "data/pdb/SCN3A_6j8e_structure_coordinates.txt"
  
  pdb_sel_gene4 = "data/pdb/SCN2A_model.pdb1"
  structure_coordinates_gene4 <- "data/pdb/SCN8A_6j8e_structure_coordinates.txt"
  
  func_colors <- c(`GABA uptake (vs wt)` = "khaki",`Total expression (vs wt)` = "indianred",`Surface expression (vs wt)` = "paleturquoise")

############### SERVER ###############
shinyServer(function(input, output, session) {
  
  # output$disclaimer <- renderValueBox({
  #   valueBox(
  #     value = tags$p("Pre-release version",style = "font-size:60%"),
  #     div("Website only online during ASHG 2021. New final version will be released during AES 2021 in December."),color = "red"
  #   )
  #   
  # })

  ############### PANEL DECISION Welcome page###############
  observeEvent(input$infoBtn, {
    updateTabsetPanel(session, "TabDisplay", selected = "infoTab")
  })
  
  observeEvent(input$familyBtn, {
    updateTabsetPanel(session, "TabDisplay", selected = "familyTab")
  })
  
  observeEvent(input$variantBtn, {
    updateTabsetPanel(session, "TabDisplay", selected = "variantTab")
  })
  
  observeEvent(input$researchBtn, {
    updateTabsetPanel(session, "TabDisplay", selected = "researchTab")
  })
  
  observeEvent(input$registryBtn, {
    updateTabsetPanel(session, "TabDisplay", selected = "registryTab")
  })
  
  ##### BASIC INFORMATION #####
  
  #can be replicated for each additional gene 
  #Factor Gene1 
  
  output$Phenotype_fac1_gene1 <- renderPlotly({
    Phenotype_fac_2.fun(basic_gene1, "Epilepsy",basic_phenotype_colors_epi)
  })
  
  output$Phenotype_fac2_gene1 <- renderPlotly({
    Phenotype_fac_2.fun(basic_gene1, "Autism",basic_phenotype_colors_autism)
  })
  
  output$Phenotype_fac3_Missense <- renderPlotly({
    Phenotype_fac_3.fun(basic_gene1, "Epilepsy_syndrome","#BEBADA","Missense")
  })
  
  output$Phenotype_fac3_PTV <- renderPlotly({
    Phenotype_fac_3.fun(basic_gene1, "Epilepsy_syndrome","#FDB462","PTV")
  })
  
  output$Phenotype_fac4_gene1 <- renderPlotly({
    Phenotype_fac_2.fun(basic_gene1, "ID_after_sz_onset",basic_phenotype_colors_ID)
  })
  
  # Numeric Gene 1
  
  output$Phenotype_num1_gene1 <- renderPlotly({
    Onset_days_1.fun(basic_gene1, basic_phenotype_num, colors_gene1)
  })
  

  output$basic_legend_plot1 <- renderPlot({
    basic_onset_legend()
  })
  
  
  ##### FOR FAMILIES #####
  
  # nothing to be calculated here # 
  
  ##### VARIANT ANALYSIS ##### 
  
  #Updates after gene change 
  
  observeEvent(input$var_gene, {

    filtered_cDNA_pos <- all_exchanges.df %>% filter(cDNA_pos == input$search_cDNA_pos, Gene == input$var_gene)
    

    updatePickerInput(session = session, inputId = "search_Allele",
                      choices = unique(filtered_cDNA_pos$Allele))

    updatePickerInput(session = session, inputId = "search_cDNA_alt",
                      choices = if(all(is.na(filtered_cDNA_pos$cDNA_alt))){NA}else{unique(filtered_cDNA_pos$cDNA_alt[!is.na(filtered_cDNA_pos$cDNA_alt)])},
                      selected = filtered_cDNA_pos$cDNA_alt[1])

    updatePickerInput(session = session, inputId = "search_AA_alt",
                      choices = all_exchanges.df %>% filter(Gene == input$var_gene, AA_pos == filtered_cDNA_pos$AA_pos[1]) %>% .$AA_alt %>% unique(),
                      selected = filtered_cDNA_pos %>% filter(cDNA_alt == cDNA_alt[1]) %>% .$AA_alt)

    updatePickerInput(session = session, inputId = "search_AA_ref",
                      choices = filtered_cDNA_pos %>% filter(cDNA_alt == cDNA_alt[1]) %>% .$AA_ref,
                      selected = filtered_cDNA_pos %>% filter(cDNA_alt == cDNA_alt[1]) %>% .$AA_ref)
  })
    
  #Updates after DNA change 
  
  observeEvent(input$search_cDNA_pos, {

    filtered_cDNA_pos <- all_exchanges.df %>% filter(cDNA_pos == input$search_cDNA_pos, Gene == input$var_gene)
    filtered_cDNA_pos_aa <- all_exchanges.df %>% filter(AA_pos == filtered_cDNA_pos$AA_pos, Gene == input$var_gene)

    updatePickerInput(session = session, inputId = "search_Allele",
                      choices = unique(filtered_cDNA_pos$Allele))

    updatePickerInput(session = session, inputId = "search_cDNA_alt",
                      choices = if(all(is.na(filtered_cDNA_pos$cDNA_alt))){NA}else{unique(filtered_cDNA_pos$cDNA_alt[!is.na(filtered_cDNA_pos$cDNA_alt)])},
                      selected = "A")

    updateNumericInputIcon(session = session, inputId = "search_AA_pos",
                      value = if(all(is.na(filtered_cDNA_pos$AA_pos))){NA}else{unique(filtered_cDNA_pos$AA_pos[!is.na(filtered_cDNA_pos$AA_pos)])})

    updatePickerInput(session = session, inputId = "search_AA_ref",
                      choices = unique(filtered_cDNA_pos_aa$AA_ref))

    updatePickerInput(session = session, inputId = "search_AA_alt",
                      choices = unique(filtered_cDNA_pos_aa$AA_alt))

  })

  observeEvent(input$search_cDNA_alt, {

    filtered_cDNA_pos <- all_exchanges.df %>% filter(cDNA_pos == input$search_cDNA_pos, Gene == input$var_gene, cDNA_alt == input$search_cDNA_alt)

    updatePickerInput(session = session, inputId = "search_AA_alt",
                      selected = unique(filtered_cDNA_pos$AA_alt))

  })
  
  #Updates after AA change
  observeEvent(input$search_AA_alt, {
    filtered_AA_alt <- all_exchanges.df %>% filter(AA_pos == input$search_AA_pos, Gene == input$var_gene, AA_alt == input$search_AA_alt)

    updatePickerInput(session = session,
                      inputId = c("get_var_type"),
                      choices = c(unique(filtered_AA_alt$Vartype)),
                      selected = c(unique(filtered_AA_alt$Vartype)))
  })

  observeEvent(input$search_AA_pos, {

    filtered_AA_pos <- all_exchanges.df %>% filter(AA_pos == input$search_AA_pos, Gene == input$var_gene)

    updatePickerInput(session = session, inputId = "search_AA_ref",
                      choices = sort(unique(filtered_AA_pos$AA_ref)))

    updatePickerInput(session = session, inputId = "search_AA_alt",
                      choices = unique(filtered_AA_pos$AA_alt))
    
  })
  
  ## Actions once the search button is pressed 
  
  varFilterInput <- reactiveValues(data=NULL)
  
  varFilterInputClinvar <- reactiveValues(data=NULL)
  
  observeEvent(input$search_var_c, {
    varFilterInput$data <- all_exchanges.df %>% filter(Gene==input$var_gene) %>% filter(cDNA_pos==input$search_cDNA_pos) %>%
      filter(Allele==input$search_Allele) %>% filter(cDNA_alt==input$search_cDNA_alt)
    
    varFilterInputClinvar$data <- clinvar.df %>% filter(Gene==input$var_gene) %>% filter(AA_pos==input$search_AA_pos) %>%
      filter(AA_ref==input$search_AA_ref) %>% filter(AA_alt==input$search_AA_alt)
    
  })
  
  observeEvent(input$search_var_p, {
    varFilterInput$data <- all_exchanges.df %>% filter(Gene==input$var_gene) %>% filter(AA_pos==input$search_AA_pos) %>%
      filter(AA_ref==input$search_AA_ref) %>% filter(AA_alt==input$search_AA_alt)
    
    varFilterInputClinvar$data <- clinvar.df %>% filter(Gene==input$var_gene) %>% filter(AA_pos==input$search_AA_pos) %>%
      filter(AA_ref==input$search_AA_ref) %>% filter(AA_alt==input$search_AA_alt)
    
  })
  
  # update cDNA input according to user search input , not 100% precise as many option may be possible 
  
   
  observeEvent(input$search_var_p, {
    updateNumericInputIcon(
      session = session,
      inputId = c("search_cDNA_pos"),
      value = c(unique(varFilterInput$data$cDNA_pos)[1])
    )

    updatePickerInput(
      session = session,
      inputId = c("search_Allele"),
      choices = c(unique(varFilterInput$data$Allele)[1])
    )

    updatePickerInput(
      session = session,
      inputId = c("search_cDNA_alt"),
      selected = c(unique(varFilterInput$data$cDNA_alt))
                   
    )
     
  })
  
  
  ##ClinVar variant interpretation
  output$ClinVarbox <- renderValueBox({
    
    clinvar_input <- ifelse(nrow(varFilterInputClinvar$data)>0,varFilterInputClinvar$data$ClinicalSignificance[1],"No ClinVar variants available")
    
    clinvar_link <-  ifelse(nrow(varFilterInputClinvar$data)>0,varFilterInputClinvar$data$Link[1],"https://www.ncbi.nlm.nih.gov/clinvar/")
    
    valueBox(
      value = tags$p(clinvar_input,style = "font-size:60%"),
      div(HTML(paste0("To access variant click <a style=color:black;  href=\"",clinvar_link,"\">here</a>"))),
      # div("To access variant click ", shiny::a("here", 
      #                                 href=clinvar_link, 
      #                                 target="_blank",
      #                                 <a style=color:white;
      #                                 )),
      color = "red"
    )
  })
  
##### Variant Information #####
  
  output$geneBox1 <- renderValueBox({
    valueBox(
      value = tags$p(paste0("Gene: ",unique(varFilterInput$data$Gene)),style = "font-size:50%"), 
      div("Transcript: ",unique(varFilterInput$data$Transcript),br(),br(),""), icon = icon("dna"),
      color = "purple"
    )
  })
  
  output$geneBox2 <- renderValueBox({
    
    if(!is.null(varFilterInput$data)){
      
      p_old <- unique(varFilterInput$data$AA_ref)  
      p_old_aa <- seqinr::a(p_old)
      p_new <- unique(varFilterInput$data$AA_alt) 
      p_new_aa <- ifelse(p_new != "Stop",seqinr::a(p_new),"*")
      p_domain <- Domain_gene.df %>% filter(Gene == varFilterInput$data$Gene[1], AA_pos == varFilterInput$data$AA_pos[1]) %>% .$Domain %>% unique()
      p_pos <- unique(varFilterInput$data$AA_pos)
      
    }else{
      
      p_old <- ""
      p_old_aa <- ""
      p_new <- ""
      p_new_aa <- ""
      p_domain <- ""
      p_pos <- ""
      
    }
    
    valueBox(
      value = tags$p(paste0("Domain: ", p_domain),style = "font-size:50%"),  
      div("Amino Acid Position: ",p_pos, br(),
          paste0("Amino Acid Change: ", p_old, " (",p_old_aa, ") "), "-" ,
          paste0(p_new, " (",ifelse(p_new != "Stop",p_new_aa,"*"),")")), icon = icon("dna"),
      color = "light-blue"
    )
  })
  
  output$geneBox3 <- renderValueBox({
    
    if(is.null(varFilterInput$data)){
      
      valueBox(
        value = tags$p(paste0("Control variants: "),style = "font-size:50%"),  
        div(paste0("gnomAD allele count at same position: "),br(),
                                                                  paste0("gnomAD allele frequency at same substitution: ")), icon = icon("dna"),
        color = "green"
      )
      
    }else{
    
      gnomad_count_exchange <- round(extract_gnomad_features(Control_data.df,varFilterInput$data,"Allele count","exchange"),2)
      
      gnomad_freq_exchange <- extract_gnomad_features(Control_data.df,varFilterInput$data,"Allele freq","exchange") 
      
      if(gnomad_freq_exchange !=0){
         gnomad_freq_exchange <- paste0(str_sub(gnomad_freq_exchange,1,3),str_sub(gnomad_freq_exchange,-3,-1))
      }
      
      gnomad_count_position <- extract_gnomad_features(Control_data.df,varFilterInput$data,"Allele count","position")
      
      valueBox(
        value = tags$p(paste0("Control variants: ",  gnomad_count_exchange),style = "font-size:50%"),  
        div(paste0("gnomAD allele count at same position: ", gnomad_count_position),br(),
                                                                  paste0("gnomAD allele frequency at same substitution ", gnomad_freq_exchange)), icon = icon("dna"),
        color = "green"
      )
    }
  })
  
  # Table with patient information ####
  output$patientTable <- DT::renderDataTable({
    
    print(varFilterInput$data)
  
    
    validate(
      need(!plyr::empty(varFilterInput$data),
           "There is no data that matches your filters.")) 
   
    datatable(Patient_data_missense_only.df %>% filter(Gene==varFilterInput$data$Gene[1]) %>%
                filter(AA_pos==varFilterInput$data$AA_pos[1]) %>%
                filter(AA_alt == varFilterInput$data$AA_alt[1],
                       !is.na(Published_in)) %>%
                select(Transcript,Gene, Domain, cDNA, Protein, Inheritance, Epilepsy, Sz_onset, Autism,Published_in),
              colnames = c("Transcript","Gene", "Domain"," cDNA", "Protein", "Inheritance",phenotype_name1,phenotype_name2,phenotype_name3,"Source"),
              options = list(dom = 't', scrollY = TRUE), escape=FALSE)
  })
  
  observe_helpers(withMathJax = TRUE)
  
  output$comparePlot <- renderPlotly({

    validate(need(
      !plyr::empty(varFilterInput$data),
      "There is no data that matches your filters."
    ))
  

    z <- Patient_data_missense_only.df %>%
      filter(case_when(
        input$compareButtons =="Variant Type" ~ Vartype ==varFilterInput$data$Vartype,
        input$compareButtons =="Protein region" ~ Domain_color==varFilterInput$data$Domain_color,
        #input$compareButtons == "Functional Consequence" ~ CALL == varFilterInput$data$CALL,
        input$compareButtons =="Amino Acid Position" ~ AA_pos==varFilterInput$data$AA_pos)) %>%
      select(Transcript,Gene, Domain, cDNA, Protein, Inheritance, Epilepsy,Sz_onset,Autism,Published_in)

    plotty1 <- plot_ly(data = z  %>%
                       dplyr::rename(selected_pheno = "Epilepsy") %>%
                       filter(selected_pheno  != "NA") %>% 
                       mutate(selected_pheno  =ifelse(selected_pheno  == "Yes"," Yes",selected_pheno )) %>% 
                       add_count(selected_pheno, name = "n_pheno") %>% distinct(selected_pheno,n_pheno) %>% 
                       assign("save",.,envir = .GlobalEnv),
                      x = ~ selected_pheno, 
                      y = ~ n_pheno, 
                      color = ~ selected_pheno, 
                      colors = basic_phenotype_colors_epi,
                      type = "bar", 
                      hoverinfo = "text", showlegend = FALSE,
                      text= ~ paste0(n_pheno, " individuals)")) %>% 
                      layout(title="Epilepsy", 
                             font=plotly_font,
                             xaxis = list(title="",showline = T, tickangle = 45),
                             yaxis = list(title="N of individuals",showline = T),
                             margin = list(b = 160)) %>%
                      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    plotty2 <- plot_ly(data = z  %>%
                    dplyr::rename(selected_pheno = "Autism") %>%
                    filter(selected_pheno  != "NA") %>% 
                    mutate(selected_pheno  =ifelse(selected_pheno  == "Yes"," Yes",selected_pheno )) %>% 
                    add_count(selected_pheno, name = "n_pheno") %>% distinct(selected_pheno,n_pheno) %>% 
                    assign("save",.,envir = .GlobalEnv),
                    x = ~ selected_pheno, 
                    y = ~ n_pheno, 
                    color = ~ selected_pheno, 
                    colors = basic_phenotype_colors_autism,
                    type = "bar", 
                    hoverinfo = "text", showlegend = FALSE,
                    text= ~ paste0(n_pheno, " individuals)")) %>% 
                    layout(title="", 
                           font=plotly_font,
                           xaxis = list(title="",showline = T, tickangle = 45),
                           yaxis = list(title="N of individuals",showline = T),
                           margin = list(b = 160)) %>%
                    config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    plotty3 <- plot_ly() %>%
                 add_boxplot(data = z,
                 x = "Seizure onset",
                 y = ~Sz_onset,
                 jitter = 0.3,
                 name = "Seizure onset",
                 boxpoints = "outliers",
                 pointpos = 0,
                 color =I("gray"),
                 marker = list (color = 'gray',
                 line = list(color ="black", width = 1)),
                 line = list(color ="black", width = 2)) %>%
                 layout(yaxis = list(type = "log",
                        title = "Seizure onset (months)",
                        tickvals = list(5,10,20,50,100,200,500,1000,2000,5000,10000),
                        tickmode = "array",
                        showline = T),
                 margin = list(b = 160),
                 xaxis = list(showline = T))
  
    subplot(plotty1,plotty2,plotty3, nrows = 1,titleY = T) %>% 
      layout(annotations = list(
        list(x = 0.1 , y = 1.1, text = "Epilepsy", showarrow = F, xref='paper', yref='paper'),
        list(x = 0.5 , y = 1.1, text = "Autism", showarrow = F, xref='paper', yref='paper'),
        list(x = 0.9 , y = 1.1, text = "Seizure onset", showarrow = F, xref='paper', yref='paper'))
      )
  })
  
  output$compare_act <- renderPlotly({
    
    var_ex.df <- Patient_data_missense_only.df %>%
      filter(AA_pos == varFilterInput$data$AA_pos, AA_alt == varFilterInput$data$AA_alt) %>% #"varFilterInput$data$AA_pos)
      filter(!is.na(`GABA uptake (vs wt)`)) %>% 
      mutate(label = "Same variant") %>% 
      distinct(label,AA_pos,AA_alt,`GABA uptake (vs wt)`)
    
    var_pos.df <- Patient_data_missense_only.df %>%
      filter(AA_pos == varFilterInput$data$AA_pos) %>% #"varFilterInput$data$AA_pos)
      filter(!is.na(`GABA uptake (vs wt)`)) %>% 
      mutate(label = "Same position")%>% 
      distinct(label,AA_pos,AA_alt,`GABA uptake (vs wt)`)
    
    var_domain.df <- Patient_data_missense_only.df %>%
      filter(Domain_color==varFilterInput$data$Domain_color) %>% #varFilterInput$data$Domain_color)
      filter(!is.na(`GABA uptake (vs wt)`)) %>% 
      mutate(label = "Same domain") %>% 
      distinct(label,AA_pos,AA_alt,`GABA uptake (vs wt)`)
    
    var_pat.df <- Patient_data_missense_only.df %>% 
      mutate(label = "All patient variants") %>% 
      distinct(label,AA_pos,AA_alt,`GABA uptake (vs wt)`)
    
    var_con.df <- Control_data.df %>% 
      left_join(functional.df %>% filter(Vartype == "Missense"), by = c("AA_pos" = "AA_pos","AA_alt" = "AA_alt","AA_ref" = "AA_ref","Vartype" = "Vartype")) %>% 
      filter(!is.na(`GABA uptake (vs wt)`)) %>% 
      mutate(label = "All control variants") %>% 
      distinct(label,AA_pos,AA_alt,`GABA uptake (vs wt)`) %>% 
      select(label,`GABA uptake (vs wt)`)
    
    
    var_all.df <- rbind(var_pos.df,var_domain.df,var_pat.df) %>% 
      select(label,`GABA uptake (vs wt)`) %>% 
      rbind(var_con.df)
    
    validate(
      need(nrow(var_all.df) >0,
           "There is no data that matches your filters.")) 
    
    plot <- plot_ly() %>%
      add_boxplot(data = var_all.df,
                  x = ~factor(label),
                  y = ~`GABA uptake (vs wt)`,
                  jitter = 0.3,
                  color = ~label,
                  boxpoints = "all",
                  pointpos = 0,
                  colors = c("black","darkred","orange","yellow"),
                  marker = list (color = 'gray',
                                 line = list(color ="black", width = 1)),
                  line = list(color ="black", width = 2)) %>% 
      add_trace(data = var_ex.df, x= ~label, y = ~`GABA uptake (vs wt)`, type = "box",
                boxpoints = "all",
                pointpos = 0,
                jitter = 0.3,
                marker = list (color = 'gray',
                               line = list(color ="purple", width = 1)),
                line = list(color ="purple", width = 2))
    
    if(var_all.df$label %>% unique() %>% length() == 5){
      plot <- plot %>% 
        layout(yaxis = list(title = "GABA uptake (vs wt)"),
               xaxis = list(title = ""),
               showlegend = FALSE,
               font=plotly_font) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    } else if(var_all.df$label %>% unique() %>% length() == 4){
      plot <- plot %>% 
        layout(yaxis = list(title = "GABA uptake (vs wt)"),
               xaxis = list(title = ""),
               margin = list(r = 100, l = 100),
               showlegend = FALSE,
               font=plotly_font) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    } else{
      plot <- plot %>% 
        layout(yaxis = list(title = "GABA uptake (vs wt)"),
               xaxis = list(title = ""),
               margin = list(r = 250, l = 250),
               showlegend = FALSE,
               font=plotly_font) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    }
    
  })
  
  output$Var_analysis_compare_var <- renderR3dmol({
    
    
    structure.df <- read_delim("data/pdb/7sk2_struc.txt",delim = "\t") %>%
      mutate(Aminoacid = aaa(Aminoacid)) %>%
      select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
    
    
    pat_var.df <- Patient_data_missense_only.df %>% 
      left_join(structure.df,by = c("AA_pos" = "Uniprot_position"))%>%
      filter(!is.na(Position_in_structure)) %>% 
      .$Position_in_structure %>% unique()
    
    con_var.df <- Control_data.df %>% 
      filter(Vartype == "Missense") %>% 
      left_join(structure.df,by = c("AA_pos" = "Uniprot_position"))%>%
      filter(!is.na(Position_in_structure)) %>% 
      .$Position_in_structure %>% unique()
    
    sel_var <- all_exchanges.df %>% 
      filter(AA_pos == varFilterInput$data$AA_pos) %>% 
      left_join(structure.df,by = c("AA_pos" = "Uniprot_position"))%>%
      filter(!is.na(Position_in_structure)) %>% 
      .$Position_in_structure %>% unique()

    sub_scale <- c(1,1)
    struc_color <- "white"
    
    rot = 270
    rot_axis = "x"
    spin_axis = "vy"
    
    
    
    modelo <- r3dmol(
      viewer_spec = m_viewer_spec(
        cartoonQuality = 10,
        lowerZoomLimit = 50,
        upperZoomLimit = 1000
      )
    )
    
    modelo <- modelo %>% m_add_model(data = "data/pdb/7sk2.pdb", format = "pdb")
    
    if(length(sel_var) != 0){
      # Zoom to encompass the whole scene
      modelo <-modelo %>% m_zoom_to() %>%
        # Set color o cartoon representation
        m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure
        # Set subunit colors
        m_set_style(
          sel = m_sel(chain = c("A")),
          style = m_style_cartoon(color = "white")
        ) %>%
        # visualize variants all
        m_set_style(
          sel = m_sel(resi = pat_var.df,
                      atom = "CA",
                      chain = c("A")),
          style = m_style_sphere(colorScheme = NULL,
                                 color = "red",
                                 scale = 0.8)
        ) %>% 
        m_set_style(
          sel = m_sel(resi = con_var.df,
                      atom = "CA",
                      chain = c("A")),
          style = m_style_sphere(colorScheme = NULL,
                                 color = "black",
                                 scale = 0.8)
        ) %>% 
        m_set_style(
          sel = m_sel(resi = sel_var,
                      atom = "CA",
                      chain = c("A")),
          style = m_style_sphere(colorScheme = NULL,
                                 color = "magenta",
                                 scale = 1.8)
        ) 
    }else{
      
      modelo <-modelo %>% m_zoom_to() %>%
        # Set color o cartoon representation
        m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure
        # Set subunit colors
        m_set_style(
          sel = m_sel(chain = c("A")),
          style = m_style_cartoon(color = "white")
        )
      
    }
    
    return(modelo)
    
  })
  
  output$Var_analyis_hotzone <- renderR3dmol({
  
    structure.df <- read_delim("data/pdb/7sk2_struc.txt",delim = "\t") %>%
      mutate(Aminoacid = aaa(Aminoacid)) %>%
      select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
    
    
    color_map.df <- read_delim("data/pdb/activity_bubble.txt",delim = "\t") %>% 
      left_join(structure.df,by = c("AA_pos" = "Uniprot_position"))%>%
      filter(!is.na(Position_in_structure)) %>% 
      distinct(Position_in_structure,act_group)
    
    sel_var <- all_exchanges.df %>% 
      filter(AA_pos %in% varFilterInput$data$AA_pos) %>% 
      left_join(structure.df,by = c("AA_pos" = "Uniprot_position"))%>%
      filter(!is.na(Position_in_structure)) %>% 
      .$Position_in_structure %>% unique()
    
    sub_scale <- c(1,1)
    struc_color <- "white"
    
    rot = 270
    rot_axis = "x"
    spin_axis = "vy"
    
    modelo <- r3dmol(
      viewer_spec = m_viewer_spec(
        cartoonQuality = 10,
        lowerZoomLimit = 50,
        upperZoomLimit = 1000
      )
    )
    
    modelo <- modelo %>% m_add_model(data = "data/pdb/7sk2.pdb", format = "pdb")
    
    if(length(sel_var) != 0){
      
      modelo <-modelo %>% m_zoom_to() %>%
        # Set color o cartoon representation
        m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure
        # Set subunit colors
        m_set_style(
          sel = m_sel(chain = c("A")),
          style = m_style_cartoon(color = "white")
        ) %>%
        # visualize variants all
        m_set_style(
          sel = m_sel(resi = color_map.df %>% filter(act_group == "<10") %>% .$Position_in_structure,
                      atom = "CA",
                      chain = c("A")),
          style = m_style_cartoon(color = "red")
        ) %>% 
        m_set_style(
          sel = m_sel(resi = color_map.df %>% filter(act_group == ">42.8") %>% .$Position_in_structure,
                      atom = "CA",
                      chain = c("A")),
          style = m_style_cartoon(color = "yellow")
        )%>% 
        m_set_style(
          sel = m_sel(resi = sel_var,
                      atom = "CA",
                      chain = c("A")),
          style = m_style_sphere(colorScheme = NULL,
                                 color = "magenta",
                                 scale = 1.8)
        ) 
      
      
    } else{
      
      modelo <-modelo %>% m_zoom_to() %>%
        # Set color o cartoon representation
        m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure
        # Set subunit colors
        m_set_style(
          sel = m_sel(chain = c("A")),
          style = m_style_cartoon(color = "white")
        ) %>%
        # visualize variants all
        m_set_style(
          sel = m_sel(resi = color_map.df %>% filter(act_group == "<10") %>% .$Position_in_structure,
                      atom = "CA",
                      chain = c("A")),
          style = m_style_cartoon(color = "red")
        ) %>% 
        m_set_style(
          sel = m_sel(resi = color_map.df %>% filter(act_group == ">42.8") %>% .$Position_in_structure,
                      atom = "CA",
                      chain = c("A")),
          style = m_style_cartoon(color = "yellow")
        )
      
    }
    
    return(modelo)
    
  })
  

  output$Var_analyis_paraz <- renderPlotly({
    
    validate(need(
      nrow(varFilterInput$data) >0,
      "There is no data that matches your filters."
    ))
    
    
    paraz_mtr_input.df <- paraz_mtr.df %>% 
      filter(Gene == varFilterInput$data$Gene[1]) %>% 
      left_join(Patient_data_missense_only.df %>% filter(Gene == varFilterInput$data$Gene[1], Vartype == "Missense") %>% mutate(p = 1) %>% select(p, AA_pos)) %>% 
      left_join(Control_data.df %>% filter(Gene == varFilterInput$data$Gene[1]) %>% mutate(g = 1) %>% select(g, AA_pos)) %>% 
      replace(is.na(.),0) %>% 
      assign("save",.,envir = .GlobalEnv) %>% 
      filter(p == 1) %>% 
      mutate(g = "pathogenic") %>% 
      rbind(save %>% filter(g == 1) %>% mutate(g = "control"))
    
    col1 <- "darkred"
    col2 <- "darkblue"
    
    plot <- plot_ly(colors=c(col1, col2)) %>% 
      add_boxplot(data = paraz_mtr_input.df %>% filter(Gene == varFilterInput$data$Gene[1], g == "pathogenic"),
                  y = ~Paraz_score, type = "box", x= 0,
                  color = I(col1)) %>% 
      add_trace(y = ~ paraz_mtr_input.df %>% filter(Gene == varFilterInput$data$Gene[1], g == "control") %>%  .$Paraz_score,
                type = "box", x = 1,
                color = I(col2)) %>% 
      add_trace(data = paraz_mtr.df%>%
                  filter(Gene == varFilterInput$data$Gene[1],
                         AA_pos == varFilterInput$data$AA_pos[1]) ,
                y = ~Paraz_score, x = 0.5,
                marker = list(color = "black",
                              size = 15),
                name = "Paraz-score",
                mode = "markers",
                type = "scatter") %>%
      add_segments(x = c(-0.5), 
                   xend = 1.5, 
                   y = 0, 
                   yend = 0,
                   line = list(color = 'gray', width = 4, dash = 'dot')) %>% 
      layout(yaxis = list(title = "Paraz-score"),
             xaxis = list(title = varFilterInput$data$Gene[1],
                          zeroline = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          showgrid = FALSE),
             font = plotly_font,
             showlegend = F) %>% 
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    
    return(plot)
    
  })
  
  
  output$paraz_legend <- renderPlot({
    
    legend <- data.frame(x=c(1,4,7), y=c(3,3,3), text=c("Pathogenic variants", "Control variants", "Selected variant"))
    plot <- ggplot(legend, aes(x=x, y=y, color=text))+
      geom_point(size = 8)+
      scale_color_manual(values = c("darkblue","darkred","black"))+
      ylim(c(2.9,3.1))+
      xlim(c(0.8,9))+
      theme_void()+
      geom_text(aes(label=text), hjust=-0.2, color="black", size =5)+
      theme(legend.position = "none")
    
    
    return(plot)
    
  })
  
  
  output$Var_analyis_mtr <- renderPlotly({
    
    validate(need(
      nrow(varFilterInput$data) >0,
      "There is no data that matches your filters."
    ))
    
    
    paraz_mtr_gene.df <- paraz_mtr.df %>% filter(Gene == varFilterInput$data$Gene[1])
    
    mtr_25 <- paraz_mtr_gene.df %>% .$MTR_score %>% sort %>% .[floor(length(.)/4)]
    mtr_5 <- paraz_mtr_gene.df %>% .$MTR_score %>% sort %>% .[floor(length(.)/20)]
    
    
    paraz_mtr_input.df <- paraz_mtr.df %>% 
      filter(Gene == varFilterInput$data$Gene[1]) %>% 
      left_join(Patient_data_missense_only.df %>% filter(Gene == varFilterInput$data$Gene[1], Vartype == "Missense") %>% mutate(p = 1) %>% select(p, AA_pos)) %>% 
      left_join(Control_data.df %>% filter(Gene == varFilterInput$data$Gene[1]) %>% mutate(g = 1) %>% select(g, AA_pos)) %>% 
      replace(is.na(.),0) %>% 
      assign("save",.,envir = .GlobalEnv) %>% 
      filter(p == 1) %>% 
      mutate(g = "pathogenic") %>% 
      rbind(save %>% filter(g == 1) %>% mutate(g = "control"))
    
    col1 <- "darkred"
    col2 <- "darkblue"
    
    plot_ly(colors=c(col1, col2)) %>% 
      add_boxplot(data = paraz_mtr_input.df %>% filter(Gene == varFilterInput$data$Gene[1], g == "pathogenic"),
                  y = ~MTR_score, type = "box", x= 0,
                  color = I(col1)) %>% 
      add_trace(y = ~ paraz_mtr_input.df %>% filter(Gene == varFilterInput$data$Gene[1], g == "control") %>%  .$MTR_score,
                type = "box", x = 1,
                color = I(col2)) %>% 
      add_trace(data = paraz_mtr_gene.df%>% 
                  filter(AA_pos == varFilterInput$data$AA_pos[1]), 
                x = 0.5,
                y = ~MTR_score, x = "",
                marker = list(color = "black",
                              size = 15),
                name = "MTR-score",
                mode = "markers",
                type = "scatter") %>%  
      add_segments(x = c(-0.5), 
                   xend = 1.5, 
                   y = mtr_25, 
                   yend = mtr_25,
                   line = list(color = '#ffb366', width = 4, dash = 'dot')) %>% 
      add_segments(x = c(-0.5), 
                   xend = 1.5, 
                   name = "Lowest 5% MTR-score",
                   y = mtr_5, 
                   yend = mtr_5,
                   line = list(color = '#e60000', width = 4, dash = 'dot')) %>% 
      add_text(y = mtr_25-0.05,
               x = 0.15,
               text = "Lowest 25% MTR-score",
               textposition = "middle right") %>% 
      add_text(y = mtr_5-0.05,
               x = 0.15,
               text = "Lowest 5% MTR-score",
               textposition = "middle right") %>% 
      layout(yaxis = list(title = "MTR-score"),
             xaxis = list(title = varFilterInput$data$Gene[1],
                          zeroline = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          showgrid = FALSE),
             font = plotly_font,
             showlegend = F) %>% 
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE) 
    
    
  })
  
  output$Var_analyis_per <- renderPlotly({
    
    validate(need(
      nrow(varFilterInput$data) >0,
      "There is no data that matches your filters."
    ))
    
    per_input.df <-per2d_for_var_analysis.df %>% 
      filter(Gene == varFilterInput$data$Gene[1]) %>% 
      left_join(Patient_data_missense_only.df %>% filter(Gene == varFilterInput$data$Gene[1], Vartype == "Missense") %>% mutate(p = 1) %>% select(p, AA_pos)) %>% 
      left_join(Control_data.df %>% filter(Gene == varFilterInput$data$Gene[1]) %>% mutate(g = 1) %>% select(g, AA_pos)) %>% 
      replace(is.na(.),0) %>% 
      assign("save",.,envir = .GlobalEnv) %>% 
      filter(p == 1) %>% 
      mutate(g = "pathogenic") %>% 
      rbind(save %>% filter(g == 1) %>% mutate(g = "control")) %>% 
      mutate(odds = log2(odds))
    
    col1 <- "darkred"
    col2 <- "darkblue"
    
    plot <- plot_ly(colors=c(col1, col2)) %>% 
      add_boxplot(data = per_input.df %>% filter(g == "pathogenic"),
                  y = ~odds, type = "box", x= 0,
                  color = I(col1)) %>% 
      add_trace(y = ~ per_input.df %>% filter(g == "control") %>%  .$odds,
                type = "box", x = 1,
                color = I(col2)) %>% 
      add_trace(data = per2d_for_var_analysis.df %>%
                  filter(AA_pos == varFilterInput$data$AA_pos[1],
                         Gene == varFilterInput$data$Gene[1]),
                y = ~odds, x = 0.5,
                marker = list(color = "black",
                              size = 15),
                name = "Fold enrichment of pathogenic variants (log2)",
                mode = "markers",
                type = "scatter") %>% 
      add_segments(x = c(-0.5), 
                   xend = 1.5, 
                   y = 0, 
                   yend = 0,
                   line = list(color = 'gray', width = 4, dash = 'dot')) %>% 
      layout(yaxis = list(title = "Fold enrichment of pathogenic variants (log2)"),
             xaxis = list(title = varFilterInput$data$Gene[1],
                          zeroline = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          showgrid = FALSE),
             font = plotly_font,
             showlegend = F) %>% 
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    
    return(plot)
    
  })
  
  #####Research #####
  # Filter for subset of variants
  res_mod <- callModule(
    module = selectizeGroupServer,
    id = "research-filters",
    data = Patient_data.df,
    vars = c("Vartype",  "AA_alt", "Domain","Epilepsy_syndrome","Autism","ID_after_sz_onset", "Hotzone_2D","PER3D")
  )
  
  res_mod_control <- callModule(
    module = selectizeGroupServer,
    id = "research-filters",
    data = Control_data.df,
    vars = c("Vartype",  "AA_alt", "Domain","Hotzone_2D","PER3D")
  )
  
  res_mod_functional <- callModule(
    module = selectizeGroupServer,
    id = "research-filters",
    data = functional_only.df,
    vars = c("Vartype",  "AA_alt", "Domain","Hotzone_2D","PER3D")
  )
  
  output$filtered_n <- renderText({
    x <- nrow(res_mod())
    x <- paste("Number of individuals:", x)
    return(x)
  })
  
  # Table with displayed variants
  output$subsetTable <- DT::renderDataTable({
    req(res_mod())

     z <- res_mod() 
     
    patient_table <- datatable(z %>% 
                                 select(Transcript,Gene, Domain, cDNA, Protein, Inheritance, Phenotype, Onset_days,functional_effect, Published_in) , 
                               extensions = "Buttons",
                               colnames = c("Transcript","Gene", "Domain"," cDNA", "Protein", "Inheritance",phenotype_name1,phenotype_name2, "Functional Consequence","Source"), 
                               escape=FALSE,
                               options = list(dom = 'Brtip',buttons = c('csv', 'excel'), pageLength=300, scrollY = "350px"))

    annotation_table <-  datatable(z %>% 
                                     select(Transcript,Gene, Domain, cDNA, Protein, Inheritance, Phenotype, Onset_days,functional_effect, Published_in), 
                                   extensions = "Buttons",
                                   colnames = c("Transcript","Gene", "Domain"," cDNA", "Protein", "Inheritance",phenotype_name1,phenotype_name2, "Functional Consequence","Source"), 
                                   escape=FALSE,
                                   options = list(dom = 'Brtip',buttons = c('csv', 'excel'), pageLength=300, scrollY = "350px"))

    if (input$patientFunSwitch == FALSE) {
      return(patient_table)
    } else {
      return(annotation_table)
    }
  })

  ### Genotype Interface

  output$Genotype_overview_plot <- renderPlotly({
    

    # 2D lolliplot with ´SLC6A1 variants
    
    research_genotype_domain.df <- all_exchanges.df %>% 
      distinct(Domain,Gene,AA_pos,Domain_color) %>% 
      dplyr::group_by(Gene,Domain,Domain_color) %>% 
      dplyr::summarise(start = min(AA_pos),
                end = max(AA_pos)) 

    g <- ggplot(data=all_exchanges.df %>% 
                  distinct(AA_pos,Gene,Domain,Domain_color) %>% 
                  left_join(res_mod() %>% 
                            filter(!is.na(AA_pos)) %>% 
                            select(Gene,Protein,AA_pos,Vartype) %>% 
                            dplyr::group_by(Protein,AA_pos,Gene,Vartype) %>% 
                            dplyr::summarise(var_count = n()) %>% 
                            ungroup() %>% 
                            dplyr::mutate(Protein_count = paste0(" ",Protein,", Variant count ",var_count)) %>% 
                            dplyr::group_by(Gene,AA_pos,Vartype) %>% 
                            dplyr::summarise(Protein_final = paste(Protein_count, collapse = ";"))) %>% 
                  dplyr::mutate(Vartype = ifelse(Vartype %in% c("Missense","PTV"),Vartype,
                                   ifelse(!is.na(Vartype),"Other",NA)))) +
          geom_segment(aes(x=AA_pos, xend=AA_pos, y=2, yend=ifelse(Vartype=="Missense", 7,8)), colour="black")+
          geom_point(aes(x=AA_pos, y=ifelse(Vartype=="Missense", 7,8), color=Vartype, text=Protein_final))+
          geom_rect(data=research_genotype_domain.df, aes(xmin=start, xmax=end, ymin=2, ymax=2.5, fill=Domain_color, text=Domain))+
          theme_classic()+
          ylim(c(1.7,10))+
          labs( x= "Amino acid sequence")+
          scale_color_manual(values = lolliplot_fill_scheme)+
          scale_fill_manual(values = lolliplot_fill_scheme)+
          #facet_grid(Gene ~ .)+
          theme(
            text = element_text(size = 10),
            axis.title.y = element_blank(),
            axis.title.x = element_blank(),
            axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.line.y = element_blank(),
            legend.position = "none"
          )
    
   
    if (input$gnomad_m == TRUE) {
      g <- g + geom_point(data=Control_data.df ,
                          size=2, color = "black", aes(x=AA_pos, y=2, alpha=0.1*Allele_count, text=paste0("Position: ",AA_pos,", Allele count: ", Allele_count)))
    }
    
    if (input$PER == TRUE) {
      
      per_sel.df <- master.df %>% 
        left_join(per2d_for_var_analysis.df) %>% 
        filter(per == "PER") %>% 
        distinct(AA_pos,Gene)
      
      g <- g + geom_point(data= per_sel.df,
                          size=1, aes(x=AA_pos, y=1.75, alpha=0.8),color = "red")
    }

    g <- ggplotly(g, tooltip = "text") %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)  %>%
      layout(title="",font=plotly_font,
             xaxis = list(title = "Amino acid sequence")
      )

  })
  
  output$Genotype_legend_plot <- renderPlot({
    
    legend <- data.frame(x=c(1,11,21,31), y=c(1, 1, 1,1), text=c("Missense", "PTV", "Control","PER"))
    plot <- ggplot(legend, aes(x=x, y=y, color=text))+
      geom_point(size = 6)+
      scale_color_manual(values = c("Missense"="#D55E00","PTV"="#0072B2","Control" ="#000000", "PER" = "red"))+
      ylim(c(0,2))+
      xlim(c(0,40))+
      theme_void()+
      geom_text(aes(label=text), hjust=-0.4, color="black")+
      theme(legend.position = "none")
    
    return(plot)
    
  })
  
  
  output$test_struc <- renderR3dmol({
    
    r3dmol(
      viewer_spec = m_viewer_spec(
        cartoonQuality = 10,
        lowerZoomLimit = 50,
        upperZoomLimit = 350
      ),
      #id = "demo",
      #elementId = "demo"
    ) %>%
      # Add model to scene
      m_add_model(data = pdb_6zsl, format = "pdb") %>%
      # Zoom to encompass the whole scene
      m_zoom_to() %>%
      # Set style of structures
      m_set_style(style = m_style_cartoon(color = "#00cc96")) %>%
      # Set style of specific selection (selecting by secondary)
      m_set_style(
        sel = m_sel(ss = "s"),
        style = m_style_cartoon(color = "#636efa", arrows = TRUE)
      ) %>%
      # Style the alpha helix
      m_set_style(
        sel = m_sel(ss = "h"), # Style alpha helix
        style = m_style_cartoon(color = "#ff7f0e")
      ) %>%
      # Rotate the scene by given angle on given axis
      m_rotate(angle = 90, axis = "y") %>%
      # Animate the scene by spinning it
      m_spin()
    
    
  })
  
  output$threeDmolGene_all_new <- renderR3dmol({
    
    
    variant.df <- res_mod() %>%
      filter(Vartype == "Missense") %>%
      mutate(label = "pathogenic") %>%
      dplyr::group_by(AA_pos,AA_ref,Gene) %>%
      dplyr::summarise(n_occ = n()) %>%
      select(AA_pos,AA_ref,n_occ,Gene)
    
    
    validate(need(
      nrow(variant.df) >0,
      "There is no data that matches your filters."
    ))

    gnomad.df <- res_mod_control() %>%
      dplyr::group_by(AA_pos,AA_ref,Gene) %>%
      dplyr::summarise(n_occ = n()) %>%
      select(AA_pos,AA_ref,n_occ,Gene)

    structure.df <- read_delim("data/pdb/7sk2_struc.txt",delim = "\t") %>%
      mutate(Aminoacid = aaa(Aminoacid)) %>%
      select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)


    variant.df <- variant.df %>%
      left_join(structure.df,by = c("AA_pos" = "Uniprot_position","AA_ref" = "Aminoacid","Gene" = "gene")) %>%
      mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not
      filter(struc_cov == "yes") %>%
      distinct(Position_in_structure,Gene) %>%
      dplyr::group_by(Position_in_structure) %>%
      dplyr::summarise(var_mut = ifelse(n() >1,"mutiple",Gene))

    gnomad.df <- gnomad.df %>%
      left_join(structure.df,by = c("AA_pos" = "Uniprot_position","AA_ref" = "Aminoacid","Gene" = "gene")) %>%
      mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>%  #this column constraints the information whether a variant can be displayed on the structure or not
      filter(struc_cov == "yes") %>%
      distinct(Position_in_structure,Gene) %>%
      dplyr::group_by(Position_in_structure) %>%
      dplyr::summarise(var_mut = ifelse(n() >1,"mutiple",Gene))


    sub_color <- c("red","black")
    sub_scale <- c(1.2,0.8)
    struc_color <- "white"

    rot = 270
    rot_axis = "x"
    spin_axis = "vy"

    #Specify yourself- color of the cartoon per subunit
    subunit_color <- c("wheat","white") #first color for GRIN1 second or GRIN2A
    # 
    # #Model for the protein complex
    
    modelo <- r3dmol(
      viewer_spec = m_viewer_spec(
        cartoonQuality = 10,
        lowerZoomLimit = 50,
        upperZoomLimit = 1000
      )
    )
    
    modelo <- modelo %>% m_add_model(data = "data/pdb/7sk2.pdb", format = "pdb")
    
    # Zoom to encompass the whole scene
    modelo <-modelo %>% m_zoom_to() %>%
      # Set color o cartoon representation
      m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure
      # Set subunit colors
      m_set_style(
        sel = m_sel(chain = c("A")),
        style = m_style_cartoon(color = subunit_color[2])
      ) %>%
      # visualize variants all
      m_set_style(
        sel = m_sel(resi = variant.df %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = sub_color[1],
                               scale = sub_scale[1])
      )



    if  (input$gnomad_m == TRUE) {

      modelo <- modelo %>% m_set_style(
        sel = m_sel(resi = gnomad.df$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = sub_color[2],
                               scale = sub_scale[2]))

    }
    return(modelo)
    
  })
  
  ##domain ernichment 
  
  output$domain_enrichment_pat_pop <- renderPlotly({
    
    
    validate(need(
      nrow(res_mod()) >0,
      "There is no data that matches your filters."
    ))
    
    
    pat_con.df <- rbind(res_mod() %>% 
                          filter(Vartype == "Missense") %>% 
                          select(AA_pos,AA_ref,AA_alt) %>% 
                          #distinct() %>% 
                          mutate(label = "Patient"),
                        res_mod_control() %>% 
                          filter(Vartype == "Missense") %>% 
                          select(AA_pos,AA_ref,AA_alt) %>% 
                          distinct() %>% 
                          mutate(label = "Control")) %>% 
      left_join(Domain_gene.df) %>% 
      group_by(Domain_color,label) %>% 
      dplyr::summarise(n = n()) %>% 
      ungroup() %>% 
      group_by(label) %>% 
      dplyr::mutate(n_out = sum(n)-n) %>% 
      ungroup() %>% 
      pivot_wider(values_from = c(n, n_out),names_from = label) %>% 
      filter(!is.na(n_Patient),
             !is.na(n_out_Patient),
             !is.na(n_Control),
             !is.na(n_out_Control))
    
    validate(need(
      nrow(pat_con.df) >1,
      "There is not sufficient data that matches your filters."
    ))
    
    
    
    domain_order <- c("N-terminal","TM1/6","Linker","TMD-other","Scaffold","EL2","EL3","EL4","C-terminal")
    domain_order <- domain_order[which(domain_order %in% pat_con.df$Domain_color)]
    
    ggplot_p <- add_odds(pat_con.df$n_Patient,pat_con.df$n_out_Patient,pat_con.df$n_Control,pat_con.df$n_out_Control,pat_con.df$Domain_color) %>% 
      mutate(sig = ifelse(pvalue>0.05,"n","y")) %>% 
      mutate(label = factor(label, levels = domain_order)) %>%   
      arrange(label) %>% 
      assign("save",.,envir = .GlobalEnv) %>% 
      ggplot(aes(y = odds, x = label,text = paste0("OR = ",round(odds ,2),"\nP-value = ",format(pvalue,digits = 2,scientific = T))))+
      geom_point(aes(color = sig), size = 3)+
      geom_errorbar(aes(ymin = lCI, ymax = uCI), size = 1, width = 0.1)+
      theme_classic(base_size = 15)+
      scale_y_log10()+
      geom_hline(aes(yintercept = 1),linetype = "dashed")+
      scale_color_manual(values = c("black","red"))+
      theme(legend.position = "none",
            panel.border = element_blank(),
            axis.text = element_text(color = "black"),
            axis.text.x = element_text(angle = 45, hjust = 1),
            axis.ticks.x = element_line(color = "black", size = 1),
            #axis.text.y = element_blank(),
            axis.ticks.y = element_blank(),
            axis.title = element_text(color = "black"),
            plot.title = element_text(color = "black", hjust = 0.5))+
      labs(x = "",
           y = "Patient variant enrichment OR",
           title = "")+
      coord_cartesian(ylim = c(min(save$lCI)*0.8,max(save$uCI)*1.2))
    
    ggplotly(ggplot_p, tooltip = "text")
    
    
  })

  ### Phenotype Interface ####

  output$research_phenotype1 <- renderPlotly({

    plot <- res_mod() %>% 
              filter(!is.na(AA_pos)) %>% 
              select(Gene,Domain) %>% 
              mutate(Gene = "Patient") %>% 
              bind_rows(Control_data.df %>% 
                        select(Domain) %>% 
                        mutate(Gene = "Control")) %>% 
              mutate(Domain = ifelse(str_detect(Domain,"Linker"),"Linker",Domain),
                     Domain = factor(Domain, levels = c("N-Terminal Cytoplasmic","Helix-1","Helix-2","Helix-3","Helix-4","Helix-5","Helix-6","Helix-7","Helix-8","Helix-9","Helix-10","Helix-11","Helix-12","Linker","C-Terminal-Cytoplasmic"))) %>% 
              dplyr::group_by(Gene, Domain) %>% 
              dplyr::summarise(n = n())  %>% 
              ungroup() %>% 
              distinct() %>% 
              plot_ly(
                x=~Domain, y=~(n), split=~Gene, type="bar", color=~Gene, alpha = 0.6, colors= Gene_colors) %>%
              layout(title=research_phenotype1_title,font=plotly_font,
                     margin = list(t =50),
                     yaxis = list(type = "log",
                                  title = "Number",
                                  ticktext = list("5", "10","20","50", "100", "200","500","1000"), 
                                  tickvals = list(1, 5,10,20,50,100,200,500,1000),
                                  showline = T,
                                  tickmode = "array"
                     ),
                     xaxis = list(title = "", showline = T, tickangle = 45)) %>% 
              config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)

    return(plot)

  })

  output$research_phenotype2 <- renderPlotly({
    plot <- res_mod() %>%  
      filter(!is.na(AA_pos)) %>% 
      dplyr::rename(Phenotype1 = "Epilepsy") %>% 
      filter(Phenotype1 != "NA") %>% 
      filter(!is.na(Phenotype1)) %>% 
      mutate(Phenotype1 = factor(Phenotype1, levels = c("Yes","No"))) %>%
      dplyr::group_by(Gene, Phenotype1) %>% dplyr::summarise(n = n())  %>% ungroup() %>% 
      plot_ly(
        x=~Phenotype1, y=~n, type="bar", alpha = 0.8, color = ~Phenotype1,colors= basic_phenotype_colors_epi) %>%
      layout(title=research_phenotype2_title ,font=plotly_font, yaxis = list(showline = T, title = "Number"),
             margin = list(t =50),
             xaxis = list(title = "",  showline = T,tickangle = 45)
      ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)

    return(plot)


  })

  output$research_phenotype3 <- renderPlotly({

    plot <- res_mod() %>%  
      dplyr::rename(Phenotype2 = Sz_onset) %>%
      select(Gene,Phenotype2) %>% 
      plot_ly(
        x=~Gene, y=~Phenotype2, type="box", color=~Gene, alpha = 0.8, colors= "purple") %>%
      layout(title=research_phenotype3_title,font=plotly_font,
             margin = list(t =50),
             xaxis = list(title = "",  tickangle = 45, showline = T),
             yaxis = list(type = "log",
                          title = "Seizure onset (months)",
                          #ticktext = list("1", "10", "100", "1000"), 
                          tickvals = list(1,5,10,20,50,100,1000,10000),
                          showline = T,
                          tickmode = "array"
             )) %>% 
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    

    return(plot)

  })
  
  output$research_phenotype4 <- renderPlotly({
    plot <- res_mod() %>%  
      filter(!is.na(AA_pos)) %>% 
      dplyr::rename(Phenotype1 = "Autism") %>% 
      filter(!is.na(Phenotype1)) %>% 
      mutate(Phenotype1 = factor(Phenotype1, levels = c("Yes","No"))) %>%
      dplyr::group_by(Gene, Phenotype1) %>% dplyr::summarise(n = n())  %>% ungroup() %>% 
      plot_ly(
        x=~Phenotype1, y=~n, type="bar", alpha = 0.8, color = ~Phenotype1,colors= basic_phenotype_colors_autism) %>%
      layout(title=research_phenotype4_title ,font=plotly_font, yaxis = list(showline = T, title = "Number"),
             margin = list(t =50),
             xaxis = list(title = "",  showline = T,tickangle = 45)
      ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    return(plot)
    
    
  })
  
  output$research_phenotype5 <- renderPlotly({
    plot <- res_mod() %>%  
      filter(!is.na(AA_pos)) %>% 
      dplyr::rename(Phenotype1 = "ID_after_sz_onset") %>% 
      dplyr::group_by(Gene, Phenotype1) %>% dplyr::summarise(n = n())  %>% ungroup() %>% 
      plot_ly(
        x=~Phenotype1, y=~n, type="bar", alpha = 0.8, color = ~Phenotype1,colors= basic_phenotype_colors_ID) %>%
      layout(title=research_phenotype5_title ,font=plotly_font, yaxis = list(showline = T, title = "Number"),
             margin = list(t =50),
             xaxis = list(title = "",  showline = T,tickangle = 45)
      ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    
    return(plot)
    
    
  })
  
  output$research_phenotype6 <- renderPlotly({
    plot <- res_mod() %>%  
      filter(!is.na(AA_pos)) %>% 
      dplyr::rename(Phenotype1 = "Epilepsy_syndrome") %>% 
      filter(Phenotype1  != "Unavailable") %>% 
      mutate(Phenotype1 = ifelse(Phenotype1 == "No seizures","No seizures",Phenotype1)) %>% 
      dplyr::group_by(Gene, Phenotype1) %>% dplyr::summarise(n = n())  %>% ungroup() %>% 
      plot_ly(
        x=~Phenotype1, y=~n, type="bar", alpha = 0.8, color = ~Phenotype1,colors= c("#bdb0b9","#5a94f4","#491919")) %>%
      layout(title=research_phenotype6_title ,font=plotly_font, yaxis = list(showline = T, title = "Number"),
             margin = list(t =50),
             xaxis = list(title = "",  showline = T,tickangle = 45)
      ) %>% config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    return(plot)
    
    
  })

  ### Research Functional Interface #####
  
  #structure 
  
  output$research_hotzone <- renderR3dmol({
    
    
    structure.df <- read_delim("data/pdb/7sk2_struc.txt",delim = "\t") %>%
      mutate(Aminoacid = aaa(Aminoacid)) %>%
      select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
    
    
    color_map.df <- read_delim("data/pdb/activity_bubble.txt",delim = "\t") %>% 
      left_join(structure.df,by = c("AA_pos" = "Uniprot_position"))%>%
      filter(!is.na(Position_in_structure)) %>% 
      distinct(Position_in_structure,act_group)
    
    
    sub_scale <- c(1,1)
    struc_color <- "white"
    
    rot = 270
    rot_axis = "x"
    spin_axis = "vy"
    
    
    
    modelo <- r3dmol(
      viewer_spec = m_viewer_spec(
        cartoonQuality = 10,
        lowerZoomLimit = 50,
        upperZoomLimit = 1000
      )
    )
    
    modelo <- modelo %>% m_add_model(data = "data/pdb/7sk2.pdb", format = "pdb")
    
    # Zoom to encompass the whole scene
    modelo <-modelo %>% m_zoom_to() %>%
      # Set color o cartoon representation
      m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure
      # Set subunit colors
      m_set_style(
        sel = m_sel(chain = c("A")),
        style = m_style_cartoon(color = "white")
      ) %>%
      # visualize variants all
      m_set_style(
        sel = m_sel(resi = color_map.df %>% filter(act_group == "<10") %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_cartoon(color = "red")
      ) %>% 
      m_set_style(
        sel = m_sel(resi = color_map.df %>% filter(act_group == ">42.8") %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_cartoon(color = "yellow")
      )
    
    return(modelo)
    
  })
  
  
  output$research_dist_struc <- renderR3dmol({
    
    structure.df <- read_delim("data/pdb/7sk2_struc.txt",delim = "\t") %>%
      mutate(Aminoacid = aaa(Aminoacid)) %>%
      select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
    
    
    color_map.df <- read_delim("data/pdb/activity_bubble.txt",delim = "\t") %>% 
      left_join(structure.df,by = c("AA_pos" = "Uniprot_position"))%>%
      filter(!is.na(Position_in_structure)) %>% 
      distinct(Position_in_structure,act_group)
    
    Tm1_res <- Domain_gene.df %>% 
      filter(Domain_color == "TM1/6",
             AA_pos < 200) %>% 
      .$AA_pos
    
    Tm6_res <- Domain_gene.df %>% 
      filter(Domain_color == "TM1/6",
             AA_pos > 200) %>% 
      .$AA_pos
    
    
    sub_scale <- c(1,1)
    struc_color <- "white"
    
    rot = 270
    rot_axis = "x"
    spin_axis = "vy"
    
    
    
    modelo <- r3dmol(
      viewer_spec = m_viewer_spec(
        cartoonQuality = 10,
        lowerZoomLimit = 50,
        upperZoomLimit = 1000
      )
    )
    
    modelo <- modelo %>% m_add_model(data = "data/pdb/GAT1_axis_pdb.pdb", format = "pdb")
    
    # Zoom to encompass the whole scene
    modelo <-modelo %>% m_zoom_to() %>%
      # Set color o cartoon representation
      m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure
      # Set subunit colors
      m_set_style(
        sel = m_sel(chain = c("A")),
        style = m_style_cartoon(color = "white")
      ) %>%
      # visualize variants all
      m_set_style(
        sel = m_sel(resi = 601,
                    chain = c("A")),
        style = m_style_stick(colorScheme = NULL,
                              color = "mediumvioletred")
      ) %>% 
      m_set_style(
        sel = m_sel(resi = 0),
        style = m_style_sphere(colorScheme = NULL,
                               scale = 0.5,
                               color = "mediumseagreen")
      ) %>% 
      m_set_style(
        sel = m_sel(resi = Tm1_res),
        style = m_style_cartoon(
          color = "lightskyblue")
      ) %>% 
      m_set_style(
        sel = m_sel(resi = Tm6_res),
        style = m_style_cartoon(
          color = "sandybrown")
      )
    modelo 
    
    return(modelo)
    
    
  })
  
  output$research_var_map1 <- renderR3dmol({
    
    if(input$pat_only == TRUE){
      
      data.df <- res_mod()
      y_label <- "Number of patients"
      titel_label <- "Functionally tested patients per domain"
      
    }else{
      
      data.df <- res_mod_functional()
      y_label <- "Number of variants"
      titel_label <- "Functionally tested variants per domain"
    }
    
    
    map_func_var(data.df,"<10","red")
    
  })
  
  output$research_var_map2 <- renderR3dmol({
    
    if(input$pat_only == TRUE){
      
      data.df <- res_mod()
      y_label <- "Number of patients"
      titel_label <- "Functionally tested patients per domain"
      
    }else{
      
      data.df <- res_mod_functional()
      y_label <- "Number of variants"
      titel_label <- "Functionally tested variants per domain"
    }
    
    
    map_func_var(data.df,"10-42.8","orange")
    
  })
  
  output$research_var_map3 <- renderR3dmol({
    
    if(input$pat_only == TRUE){
      
      data.df <- res_mod()
      y_label <- "Number of patients"
      titel_label <- "Functionally tested patients per domain"
      
    }else{
      
      data.df <- res_mod_functional()
      y_label <- "Number of variants"
      titel_label <- "Functionally tested variants per domain"
    }
    
    
    map_func_var(data.df,">42.8","yellow")
    
  })

  output$research_functional1 <- renderPlotly({
    
    
    if(input$pat_only == TRUE){
      
      data.df <- res_mod()
      y_label <- "Number of patients"
      titel_label <- "Functionally tested patients per domain"
      
    }else{
      
      data.df <- res_mod_functional()
      y_label <- "Number of variants"
      titel_label <- "Functionally tested variants per domain"
    }
    
    validate(
      need(nrow(data.df)>0,
           "There is no data that matches your filters.")) 
    
    data.df %>% 
      filter(Vartype == "Missense") %>% #colnames()
      group_by(Domain_color) %>% 
      dplyr::summarise(`GABA uptake (vs wt)` = sum(!is.na(`GABA uptake (vs wt)`)),
                       `Surface expression (vs wt)`  = sum(!is.na(`Surface expression (vs wt)`)),
                       `Total expression (vs wt)` = sum(!is.na(`Total expression (vs wt)`))) %>% 
      mutate(Domain_color = factor(Domain_color, levels = c("N-terminal","TM1/6","Linker","TMD-other","Scaffold","EL2","EL3","EL4","C-terminal"))) %>% 
      ungroup() %>% 
      pivot_longer(cols = c(`GABA uptake (vs wt)`, `Surface expression (vs wt)`, `Total expression (vs wt)`)) %>% 
      plot_ly(
        x=~Domain_color , y=~(value), split=~name, type="bar", color=~name, alpha = 0.6, colors= func_colors) %>%
      layout(title=titel_label,font=plotly_font,
             margin = list(t =50),
             yaxis = list(type = "log",
                          title = y_label,
                          ticktext = list("5", "10","20","50", "100", "200","500","1000"), 
                          tickvals = list(1, 5,10,20,50,100,200,500,1000),
                          showline = T,
                          tickmode = "array"
             ),
             xaxis = list(title = "", showline = T, tickangle = 45)) %>% 
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
  })
  
  
  output$research_functional2 <- renderPlotly({
    
    if(input$pat_only == TRUE){
      
      data.df <- res_mod()
      y_label <- "Number of patients"
      titel_label <- "Functionally tested patients per domain"
      
    }else{
      
      data.df <- res_mod_functional()
      y_label <- "Number of variants"
      titel_label <- "Functionally tested variants per domain"
    }
    
    validate(
      need(nrow(data.df)>0,
           "There is no data that matches your filters.")) 
    
    plot_ly(data.df %>% 
              mutate(Domain_color = factor(Domain_color, levels = c("N-terminal","TM1/6","Linker","TMD-other","Scaffold","EL2","EL3","EL4","C-terminal"))),
            x = ~Domain_color, y = ~`GABA uptake (vs wt)`,
            type = "box",
            boxpoints = "all",
            jitter = 0.3,
            pointpos = 0, 
            marker = list(color = "black"),
            line = list(color = "black"),
            fillcolor = "khaki") %>% 
      layout(title="GABA uptake per domain", 
             font=plotly_font,
             xaxis = list(title="", tickangle = 45),
             margin = list(t = 50)) %>%
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    
  })
  
  output$research_functional3 <- renderPlotly({
    
    
    if(input$pat_only == TRUE){
      
      data.df <- res_mod()
      y_label <- "Number of patients"
      titel_label <- "Functionally tested patients per domain"
      
    }else{
      
      data.df <- res_mod_functional()
      y_label <- "Number of variants"
      titel_label <- "Functionally tested variants per domain"
    }
    
    func_dist(data.df,"dist_TGI","Distance from Tiagabine vs. GABA uptake rate","Distance from ligand Tiabgabine","mediumvioletred")
    
    
  })
  
  output$research_functional4 <- renderPlotly({
    
    
    if(input$pat_only == TRUE){
      
      data.df <- res_mod()
      y_label <- "Number of patients"
      titel_label <- "Functionally tested patients per domain"
      
    }else{
      
      data.df <- res_mod_functional()
      y_label <- "Number of variants"
      titel_label <- "Functionally tested variants per domain"
    }
    
    func_dist(data.df,"dist_TM6","Distance from TM6 vs. GABA uptake rate","Distance from TM6","sandybrown")
    
    
  })
  
  output$research_functional5 <- renderPlotly({
    
    
    if(input$pat_only == TRUE){
      
      data.df <- res_mod()
      y_label <- "Number of patients"
      titel_label <- "Functionally tested patients per domain"
      
    }else{
      
      data.df <- res_mod_functional()
      y_label <- "Number of variants"
      titel_label <- "Functionally tested variants per domain"
    }
    
    func_dist(data.df,"dist_TM1","Distance from TM1 vs. GABA uptake rate","Distance from TM1","lightskyblue")
    
    
  })
  
  output$research_functional6 <- renderPlotly({
    
    
    if(input$pat_only == TRUE){
      
      data.df <- res_mod()
      y_label <- "Number of patients"
      titel_label <- "Functionally tested patients per domain"
      
    }else{
      
      data.df <- res_mod_functional()
      y_label <- "Number of variants"
      titel_label <- "Functionally tested variants per domain"
    }
    
    func_dist(data.df,"GAT1_axis_dist","Distance from GAT1-axis vs. GABA uptake rate","Distance from GAT1-axis","mediumseagreen")
    
    
  })
  
  output$research_functional7 <- renderPlotly({
    
    if(input$pat_only == TRUE){
      
      data.df <- res_mod()
      y_label <- "Number of patients"
      titel_label <- "Functionally tested patients per domain"
      
    }else{
      
      data.df <- res_mod_functional()
      y_label <- "Number of variants"
      titel_label <- "Functionally tested variants per domain"
    }
    
    data.df <- data.df %>% 
      dplyr::rename(func_effect = "Surface expression (vs wt)") %>% 
      filter(!is.na(func_effect)) %>% 
      filter(Vartype == "Missense") %>% 
      distinct(func_effect,AA_pos,AA_alt)
    
    validate(need(
      nrow(data.df) >0,
      "There is no data that matches your filters."
    ))
    
    
    plot <- plot_ly(data = data.df ,
                    y = ~func_effect, type = "box", x= 0, boxpoints = "all", pointpos = 0,
                    fillcolor = "paleturquoise",
                    marker = list(color = "black"),
                    line = list(color = "black")) %>% 
      layout(yaxis = list(title = "Surface expression (vs wt)"),
             title = "Surface expression",
             titlefont = list(size = 20),
             xaxis = list(title = "",
                          zeroline = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          showgrid = FALSE),
             font = plotly_font,
             showlegend = F,
             margin = list(t = 60)) %>% 
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    
    return(plot)
  })
  
  output$research_functional8 <- renderPlotly({
    
    if(input$pat_only == TRUE){
      
      data.df <- res_mod()
      y_label <- "Number of patients"
      titel_label <- "Functionally tested patients per domain"
      
    }else{
      
      data.df <- res_mod_functional()
      y_label <- "Number of variants"
      titel_label <- "Functionally tested variants per domain"
    }
    
    data.df <- data.df %>% 
      dplyr::rename(func_effect = "Total expression (vs wt)") %>% 
      filter(!is.na(func_effect)) %>% 
      filter(Vartype == "Missense") %>% 
      distinct(func_effect,AA_pos,AA_alt)
    
    validate(need(
      nrow(data.df) >0,
      "There is no data that matches your filters."
    ))
    
    
    plot <- plot_ly() %>% 
      add_boxplot(data = data.df ,
                  y = ~func_effect, type = "box", x= 0, boxpoints = "all", pointpos = 0,
                  fillcolor = "indianred",
                  marker = list(color = "black"),
                  line = list(color = "black")) %>% 
      layout(yaxis = list(title = "Total expression (vs wt)"),
             title = "Total expression",
             titlefont = list(size = 20),
             xaxis = list(title = "",
                          zeroline = FALSE,
                          showline = FALSE,
                          showticklabels = FALSE,
                          showgrid = FALSE),
             font = plotly_font,
             showlegend = F,
             margin = list(t = 60)) %>% 
      config(modeBarButtonsToRemove = goodbye, displaylogo = FALSE)
    
    
    return(plot)
  })
  
  
  output$threeDmolfunctional <- renderR3dmol({
    
    
    selection_data.df <- res_mod() 
    
    if(selection_data.df$Gene %>% unique() %>% length()== 4){
      
      genes_selected <-master.df$Gene %>% unique()
    }else{
      genes_selected <- selection_data.df$Gene %>% unique()
    }
    
    
    variant.df <- master.df %>% 
      left_join(Functional_data_mermer.df %>% distinct(Gene,AA_pos,functional_effect), by = c("AA_pos" = "AA_pos","Gene" = "Gene")) %>% 
      left_join(all_exchanges.df %>% distinct(Domain,Gene,AA_pos), by = c("AA_pos" = "AA_pos","Gene" = "Gene")) %>% 
      mutate(AA_ref = aaa(AA_ref)) %>% 
      filter(Domain %in% selection_data.df$Domain,
             AA_ref %in% selection_data.df$AA_ref[!is.na(selection_data.df$AA_ref)],
             Gene %in% genes_selected,
             functional_effect %in% selection_data.df$functional_effect) %>% 
      filter(!is.na(functional_effect)) %>% 
      distinct(Aln_pos,Domain,functional_effect) %>% 
      dplyr::group_by(Aln_pos) %>% 
      dplyr::summarise(functional_effect = ifelse(unique(functional_effect) %>% length() ==1,functional_effect,"complex")) %>% 
      left_join(master.df %>% filter(Gene == "SCN2A") %>% distinct(AA_pos,AA_ref,Aln_pos,Gene)) 
    
    
    structure.df <- read_delim("data/pdb/6j8e_structure_coordinates.txt",delim = "\t") %>%
      select(Uniprot_position,Aminoacid,Position_in_structure,gene,chain)
    
    
    variant.df <- variant.df %>%
      left_join(structure.df,by = c("AA_pos" = "Uniprot_position","AA_ref" = "Aminoacid","Gene" = "gene")) %>%
      mutate(struc_cov = ifelse(is.na(Position_in_structure),"no","yes")) %>% 
      filter(struc_cov == "yes") 
    
    
    sub_color <- c("#cc0000","#66ffff", "#ff8000","#6600cc","#c0c0c0")
    sub_scale <- c(1.2,0.8)
    struc_color <- "white"
    
    rot = 270
    rot_axis = "x"
    spin_axis = "vy"
    
    #Specify yourself- color of the cartoon per subunit
    subunit_color <- c("wheat","white") #first color for GRIN1 second or GRIN2A
    
    #Model for the protein complex
    
    modelo <- r3dmol(
      viewer_spec = m_viewer_spec(
        cartoonQuality = 10,
        lowerZoomLimit = 50,
        upperZoomLimit = 1000
      )
    )
    
    modelo <- modelo %>% m_add_model(data = "data/pdb/SCN2A_model.pdb1", format = "pdb")
    
    # Zoom to encompass the whole scene
    modelo <-modelo %>% m_zoom_to() %>%
      # Set color o cartoon representation
      m_set_style(style = m_style_cartoon(color = struc_color)) %>% # select a color of the structure
      # Set subunit colors
      m_set_style(
        sel = m_sel(chain = c("A")),
        style = m_style_cartoon(color = subunit_color[2])
      ) %>%
      # visualize variants all
      m_set_style(
        sel = m_sel(resi = variant.df %>% filter(functional_effect == "LoF") %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = sub_color[1],
                               scale = sub_scale[1])
      ) %>% 
      #visualize variants SCN1A
      m_set_style(
        sel = m_sel(resi = variant.df %>% filter(functional_effect == "GoF") %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = sub_color[2],
                               scale = sub_scale[1])
      ) %>% 
      #visualize variants SCN2A
      m_set_style(
        sel = m_sel(resi = variant.df %>% filter(functional_effect == "Mixed") %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = sub_color[3],
                               scale = sub_scale[1])
      )  %>% 
      # visualize variants SCN3A
      m_set_style(
        sel = m_sel(resi = variant.df %>% filter(functional_effect == "complex") %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = sub_color[4],
                               scale = sub_scale[1])
      ) %>% 
      # visualize variants SCN8A
      m_set_style(
        sel = m_sel(resi = variant.df %>% filter(functional_effect == "STW") %>% .$Position_in_structure,
                    atom = "CA",
                    chain = c("A")),
        style = m_style_sphere(colorScheme = NULL,
                               color = sub_color[5],
                               scale = sub_scale[1])
      ) 
    
    
    return(modelo)
    
  })
  
})
