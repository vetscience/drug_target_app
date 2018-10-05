library(tidyverse)
library(rlang)

# Calculation section ----
shdf <- as.tibble(read.table("merged_Sh_ranking_final.tsv", header = T, sep = "\t")) %>% rename(sp_simil=sp_sim)
shdf$MS3_id <- as.character(shdf$MS3_id)
numprot <- shdf %>% select(MS3_id) %>% count()

chembl <- as.tibble(read.table("ChEMBL_Sh_June18.tsv", header = T, sep = "\t", comment.char = "", quote = ""))
chembl$pat_exp <- as.Date(chembl$pat_exp, format = "%Y-%m-%d")
chembl$first_approv <- as.Date(chembl$first_approv, format = "%Y")
chembl <- mutate_if(chembl, is.factor, funs(as.character(.)))
chembl <- chembl %>% mutate(ro5_viol = str_replace(ro5_viol, "^None$" , NA_character_)) %>%
  mutate(ro3_pass = str_replace(ro3_pass, "^None$" , NA_character_)) %>%
  mutate(cmpd_syn = str_replace(cmpd_syn, "^None$" , NA_character_)) %>%
  mutate(nat_prod = case_when(
    nat_prod == -1 ~ "no",
    nat_prod == 1 ~ "yes",
    nat_prod == 0 ~ "no",
    TRUE ~ as.character(nat_prod)
  )
  )

chembl$ro5_viol <- as.integer(chembl$ro5_viol)

#Read in screening data for Schistosoma from ChEMBL
#Then transform every set of link_name and link into a HTML link tag, then put all links into one field for each ChEMBL value
screenS <- read_tsv(file = "schisto_screen_ChEMBL_refs.tsv", col_names = T, col_types = "ccc", comment = "") %>%
  mutate(ref_link = paste('<a href="',link,'" target="_blank" class="btn btn-primary">', link_name, '</a>', sep = "")) %>%
  select(chembl_id, ref_link) %>%  distinct() %>% group_by(chembl_id) %>% mutate(link = paste0(ref_link, collapse = " ")) %>% ungroup() %>%
  select(chembl_id, link) %>% rename(chembl_cmpd_id=chembl_id)

drugb <- as.tibble(read.table("DrugBank_Sh_final.tsv", header = T, sep = "\t", comment.char = "", quote = ""))
drugb <- mutate_if(drugb, is.factor, funs(as.character(.)))

wght <- rep(1, ncol(shdf)+2)
names(wght) <- c(colnames(shdf), "one_cmpd", "more_cmpds")
more_cmpds <- 5
ipro <- read_tsv(file = "IPRO_IDs.tsv", col_names = T, col_types = "ccc", comment = "")
prot_types <- ipro %>% select(., ipro_id, ipro_descr) %>% filter(ipro_id != "") %>% distinct()
ipro <- shdf %>% select(MS3_id) %>% left_join(.,ipro, by = "MS3_id") %>% replace_na(list(ipro_id = "", ipro_descr = "")) %>% distinct()
ipro <- ipro %>% count(ipro_id) %>% filter(n == 1) %>% left_join(ipro, ., by ="ipro_id") %>% replace_na(list(n=0)) %>% rename(uniqIpro = n)
ptdd <- c("All identifiers",paste(prot_types$ipro_id, prot_types$ipro_descr, sep = " "))
equalWeight <- 5

#set preselected cut-offs for blast coverage and similarity to determine homology
simil_cutoff <- 80
cov_cutoff <- 50

cutOff <- function(df, cut=0)
{
  if_else(df >= cut, 1, 0)
}


chembl <- chembl[chembl$target_type == "SINGLE PROTEIN",]
chembl$MS3_id <- as.character(chembl$MS3_id)


getScore <- function(x, y, c, d) {
  
  make_empty <- 0
  init_if <- 0
  homo_loop <- 0
  un11 <- 0
  un12 <- 0
  un13 <- 0
  un14 <- 0
  unique2 <- 0
  last_loop <- 0
  bind <- 0
  
  #create empty data frame with scoring columns for each criterion
  male_tpm <- rep(0,nrow(x))
  chokepoint <- rep(0,nrow(x))
  mus_leth <- rep(0,nrow(x))
  dmel_leth <- rep(0,nrow(x))
  cel_leth <- rep(0,nrow(x))
  hsa_lower_75_qtl <- rep(0,nrow(x))
  sp_cov <- x$sp_cov
  sman_cov <- x$sman_cov
  sjap_cov <- x$sjap_cov
  opvi_cov <- x$opvi_cov
  closi_cov <- x$closi_cov
  fahe_cov <- x$fahe_cov
  cel_cov <- x$cel_cov
  dmel_cov <- x$dmel_cov
  mus_cov <- x$mus_cov
  sp_simil <- x$sp_simil
  sman_simil <- x$sman_simil
  sjap_simil <- x$sjap_simil
  opvi_simil <- x$opvi_simil
  closi_simil <- x$closi_simil
  fahe_simil <- x$fahe_simil
  cel_simil <-  x$cel_simil
  dmel_simil <- x$dmel_simil
  mus_simil <- x$mus_simil
  ChCmpdOne <- rep(0,nrow(x))
  ChCmpdMore <- rep(0,nrow(x))
  DBCmpdOne <- rep(0,nrow(x))
  DBCmpdMore <- rep(0,nrow(x))
  scorenames <- x[,"MS3_id"]
  
  scores <- tibble(
    male_tpm,
    chokepoint,
    mus_leth,
    dmel_leth,
    cel_leth,
    hsa_lower_75_qtl,
    sp_cov,
    sman_cov,
    sjap_cov,
    opvi_cov,
    closi_cov,
    fahe_cov,
    cel_cov,
    dmel_cov,
    mus_cov,
    sp_simil,
    sman_simil,
    sjap_simil,
    opvi_simil,
    closi_simil,
    fahe_simil,
    cel_simil,
    dmel_simil,
    mus_simil,
    ChCmpdOne,
    ChCmpdMore,
    DBCmpdOne,
    DBCmpdMore
  )
  
  scores <- bind_cols(scorenames, scores)
  
  for (prot in c(1:nrow(x))) {
    cat(prot, "\n")
    if(x[prot, "male_tpm"] > 0 || x[prot, "female_tpm"] > 0)
    {
      scores[prot,"male_tpm"] <- y["male_tpm"]
    }
    
    #chokepoint
    if(x[prot, "chokepoint"] == TRUE){scores[prot, "chokepoint"] <- y["chokepoint"]}
    
    #lethal mus
    if(x[prot, "mus_leth"] != ""){scores[prot, "mus_leth"] <- y["mus_leth"]}
    
    #lethal dmel
    if(x[prot, "dmel_leth"] != ""){scores[prot, "dmel_leth"] <- y["dmel_leth"]}
    
    #lethal cel
    if(x[prot, "cel_leth"] != ""){scores[prot, "cel_leth"] <- y["cel_leth"]}
    
    #hsa_lower_75_qtl
    if(is.na(x[prot, "hsa_lower_75_qtl"]) || x[prot, "hsa_lower_75_qtl"] == TRUE){scores[prot, "hsa_lower_75_qtl"] <- y["hsa_lower_75_qtl"]}
    
    chembl_cmpds <- c %>% filter(MS3_id == as.character(x[prot, 1]))
    numcmpds <- chembl_cmpds %>% summarise(n = n_distinct(chembl_cmpd_id)) %>% pull()
    
    if(numcmpds > 0){scores[prot, "ChCmpdOne"] <- y["one_cmpd"]}
    
    if(numcmpds >= more_cmpds){scores[prot, "ChCmpdMore"] <- y["more_cmpds"]}
    #NTH: add number of compounds as column so one can filter on it in the output data
    
    drugb_cmpds <- d %>% filter(MS3_id == as.character(x[prot, 1]))
    numcmpds <- nrow(unique(drugb_cmpds[,"cmpd_id"]))
    
    
    if(numcmpds > 0)
    {
      scores[prot, "DBCmpdOne"] <- y["one_cmpd"]
    }
    
    if(numcmpds >= more_cmpds)
    {
      scores[prot, "DBCmpdMore"] <- y["more_cmpds"]
    }
    
    
  }
  
  scores <- select(ipro, MS3_id, uniqIpro) %>%
    group_by(MS3_id) %>%
    summarise(uniqIpro = max(uniqIpro)) %>%
    left_join(scores, ., by = "MS3_id") %>%
    replace_na(., list(uniqIpro=0)) %>% ungroup()
  
  return(scores)
}

result <- getScore(shdf, wght, chembl, drugb)
write_delim(result, path = "result.tsv", delim = "\t")
#write_delim(chembl, path = "chembl3.tsv", delim = "\t")
#write_delim(drugb, path = "drugb3.tsv", delim = "\t")
#write_delim(ipro, path = "ipro3.tsv", delim = "\t")
#write_delim(prot_types, path = "prot_types3.tsv", delim = "\t")
#write_delim(screenS, path = "screenS3.tsv", delim = "\t")
#write_delim(shdf, path = "shdf3.tsv", delim = "\t")
