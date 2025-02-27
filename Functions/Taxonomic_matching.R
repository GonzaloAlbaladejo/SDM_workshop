##----////\/</\>/>\><>/\><><>>\/><></\/\/\/\><><><|><|><>||||~/\/\/\/\/\/\#
## Extract the taxonomic information and list of sinonyms for a given specie
##----////\/</\>/>\><>/\><><>>\/><></\/\/\/\--><><><><><|>|<|<|><|>|>|>|>|~#
#

retrieve_syns<-function(spp_name,   # [Character] The species name from which to collect taxonomic information
                        n_times=100,  # [Numeric] Number of times the search is repeated until a data is found,default value = 1
                        Gbif=FALSE  # [Logical] Should we check Gbif for a taxonomic macthing of the species
)
{
  # 0. Load the packages
  list.of.packages<-c("tidyr","rredlist","taxize","data.table","stringr",
                      "rgbif","data.table","dplyr")
  
  new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
  if(length(new.packages)) install.packages(new.packages)
  
  lapply(list.of.packages,require,character.only=TRUE)
  rm(list.of.packages,new.packages)
  
  # a. Check the species names, if the name is binomial it runs the query to collect the taxonomic information
  # a.1 Removes special characters and fix the capitalization of the name
  spp.x <- stringr::str_trim(spp_name) %>% gsub(pattern = "[[:punct:]]", replacement = " ") %>% stringr::str_trim("both") %>% gsub(pattern = "  ",replacement = "") # Remove white spaces at the start and end of the string
  
  # Resolve capitalization
  CapSp <- function(x) {
    s <- strsplit(x, " ") %>% unlist()
    
    if(length(s)>1){
      paste(paste0(toupper(substring(s[1], 1,1)), substring(s[1], 2)),
            paste(tolower(s[-1]),collapse = " "), sep = " ")
      
    }else{
      paste0(toupper(substring(s[1], 1,1)), substring(s[1], 2))
    }
  }
  
  spp.x<-CapSp(spp.x)
  
  # a.2 Check if the name is related with a species/sub-specie or other taxonomic class (Class, Order, Family) 
  #     by analyzing the number of terms of the character string  
  correct_name <- NULL
  t_11 <- 1
  
  while(is.null(correct_name) && t_11 <= n_times) {
    
    try(correct_name<-gnr_resolve(sci =spp.x,
                                  resolve_once = FALSE,
                                  canonical=TRUE,
                                  best_match_only = TRUE,
                                  highestscore = TRUE),
        silent = TRUE)
    t_11 <- t_11 + 1
  }
  rm(t_11)  
  
  # Summarize the results of the name checking and correction
  if (!is.null(correct_name)){
    if(nrow(correct_name)!=0){
      y.d<-cbind(spp.x,correct_name[,colnames(correct_name)%in%
                                      c("data_source_title","score","matched_name2")]) #
      
      names(y.d)[1]<-"or_name"
      
      spp.x<-y.d$matched_name2  # Use the corrected name for the rest of taxnomic querys
      
      ifelse(y.d$matched_name2==y.d$or_name,
             y.d$Status<-"Correct",
             y.d$Status<-"Incorrect")
      
    }} else {
      y.d<-data.frame(or_name=spp.x,
                      matched_name2=NA,
                      Status="Not_found",
                      data_source_title=NA,
                      score=NA)
    }
  
  # b.Use the corrected or original name to look for taxonomic data----
  #   b.1. Get the basic data from the IUCN red list----
  #
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  sp_long <- strsplit(spp.x,split=" ") %>% unlist()
  
  g <- NULL
  t_11 <- 1
  
  while( is.null(g) && t_11 <= n_times){
    
    # try(g <- rl_search(name = spp.x,silent=TRUE)$result)
    try(g <- rl_species(genus=sp_long[1],species=sp_long[2],silent=TRUE),silent=T)
    t_11 <- t_11 + 1
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # b <- NULL
  # t_2 <- 1
  # while( is.null(b) && t_2 <= n_times) {
  #   
  #   try(b <- rl_synonyms(name=spp.x,silent=TRUE))
  #   t_2 <- t_2 + 1
  # }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
 # if (length(g)==0 & b$count==0) {
   if (is.null(g)) {  
    IUCN_id  <- NA
    IUCN_name<- NA
    IUCN_Phylum   <- NA
    IUCN_Class    <- NA
    IUCN_Order    <- NA
    IUCN_Family   <- NA
    
    IUCN_Category <- NA
    IUCN_Present <- "No"
    
    IUCN_N_syn <- NA
    IUCN_syn <- NA
    
  } else {
    # if the name of the species correspond to a synonim in the IUCN-red list, 
    # use the accepted name to retrieve the species information
    
    # if(b$count!=0 & length(g)==0){ 
    #   
    #   g <- NULL
    #   t_11 <- 1
    #   
    #   while( is.null(g) && t_11 <= n_times){
    #     
    #     try(g <- rl_search(name = b$result$accepted_name[1],silent=TRUE)$result)
    #     t_11 <- t_11 + 1
    #   }
    # }
    
    # IUCN_id  <- unique(g$taxonid) 
    IUCN_id  <- unique(g[["taxon"]][["sis_id"]])
    
    # IUCN_name <- unique(g$scientific_name)
    IUCN_name <- unique(g[["taxon"]][["scientific_name"]])
    
    # if (is.null(g$phylum)){ IUCN_Phylum <- NA } else {IUCN_Phylum <- unique(g$phylum) }
    if (length(g[["taxon"]][["phylum_name"]])==0){ IUCN_Phylum <- NA } else {IUCN_Phylum <- unique(g[["taxon"]][["phylum_name"]]) }
    
    # if (is.null(g$class)){ IUCN_Class <- NA }else {IUCN_Class <- unique(g$class) }
    if (length(g[["taxon"]][["class_name"]])==0){ IUCN_Class <- NA }else {IUCN_Class <- unique(g[["taxon"]][["class_name"]]) }
    
    # if (is.null(g$order)){ IUCN_Order <- NA }else {IUCN_Order <- unique(g$order) }
    if (length(g[["taxon"]][["order_name"]])==0){ IUCN_Order <- NA }else {IUCN_Order <- unique(g[["taxon"]][["order_name"]]) }
    
    # if (is.null(g$family)){ IUCN_Family <- NA } else {IUCN_Family <- unique(g$family)}
    if (length(g[["taxon"]][["family_name"]])==0){ IUCN_Family <- NA } else {IUCN_Family <- unique(g[["taxon"]][["family_name"]])}
    
    id_assesment <- g[["assessments"]][["assessment_id"]]
    
    if(length(id_assesment)!=0){
      a <- rl_assessment(id=id_assesment[length(id_assesment)])
      IUCN_Category <- a[["red_list_category"]][["code"]]
      
          }else{
      IUCN_Category <- NA
      }
    
    IUCN_Present <-"Yes"
    
    if(length(g[["taxon"]][["synonyms"]])!=0){
      IUCN_syn <- paste(paste(g[["taxon"]][["synonyms"]]$genus_name,g[["taxon"]][["synonyms"]]$species_name,sep=" "),collapse=";")
      IUCN_N_syn <- nrow(g[["taxon"]][["synonyms"]]) # number of sp sinonims, Subsp excluded
      
      } else {
      IUCN_N_syn<-NA
      IUCN_syn<-NA
    }
  }
  
  # b.2 Combine the IUCN information----
  IUCN_data<-data.frame(IUCN_Present,IUCN_id,IUCN_name,
                        IUCN_Category,IUCN_N_syn,IUCN_syn,
                        IUCN_Phylum,IUCN_Class,
                        IUCN_Order,IUCN_Family)
  
  # b.3. Get the basic data from the ITIS red list----
  # b.3.1 Get TSN a reference number that we are going to need to gather information from ITIS----
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  TSN <- NULL
  t_4 <- 1
  
  while( is.null(TSN) && t_4 <= n_times){
    try(TSN <- get_tsn_(spp.x,searchtype = "scientific",silent=TRUE))
    t_4 <- t_4 + 1
  }
  #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  
  if (is.null(TSN[[1]])) { 
    
    ITIS_Present<-"No"
    ITIS_id<-NA
    ITIS_name<-NA
    ITIS_is_valid<-NA
    ITIS_Phylum<-NA
    ITIS_Class<-NA
    ITIS_Order<-NA
    ITIS_Family<-NA
    
    ITIS_accept<-NA
    ITIS_syn<-NA
    ITIS_N_syn<-NA
    
    
  } else {
    
    tsn<-TSN[[1]]
    
    tsn_n<-tsn[tsn$scientificName==spp.x & tsn$nameUsage=="valid",]
    
    if(length(tsn_n$tsn)==0){ # if the original species name is not present in the returned list of synonims we igore the results
      
      ITIS_Present<-"Wrong"
      ITIS_id<-NA
      ITIS_name<-NA
      ITIS_is_valid<-NA
      ITIS_Phylum<-NA
      ITIS_Class<-NA
      ITIS_Order<-NA
      ITIS_Family<-NA
      
      ITIS_accept<-NA
      ITIS_syn<-NA
      ITIS_N_syn<-NA
      
    }else{
      
      ITIS_Present<-"Yes"
      
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      AceptName <- NULL
      t_4 <- 1
      
      while( is.null(AceptName) && t_4 <= n_times){
        
        try(AceptName <- itis_acceptname(tsn_n$tsn,silent=TRUE))
        t_4 <- t_4 + 1
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      
      if (AceptName$submittedtsn==AceptName$acceptedtsn){ # is the submitted Or_name name accepted by ITIS
        ITIS_name<-spp.x
        ITIS_accept<-"yes"
        
      }else{
        ITIS_name<-"No"                        
        ITIS_accept<-"No"
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      Syn <- NULL
      t_5 <- 1
      
      while( is.null(Syn) && t_5 <= n_times){
        
        try(Syn <- synonyms(tsn$tsn,db="itis"))
        t_5 <- t_5 + 1
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      
      # Syn_l is a list of possible codes with their respective names
      # We need to unify the list and extract the names
      
      Syn<-rbindlist(Syn,fill=TRUE)
      # count the number of FALSE and TRUE       #strsplit(, " "), length
      z<-as.character(is.na(Syn$syn_name))
      # if the length equals cero, there is no syn names
      
      if (length(z[z==TRUE])==0){
        ITIS_N_syn<-NA
        ITIS_syn<-NA
      } else {
        
        z_b<-as.character(sapply(strsplit(Syn$syn_name,split = " "),length)==2)
        
        ITIS_N_syn<-length(z_b[z_b==TRUE]) # number of sp sinonims, Subsp excluded
        ITIS_syn<-paste(Syn$syn_name[sapply(strsplit(Syn$syn_name,split = " "),length)==2],
                        collapse=";") # combine the names into a single string with all the sinonims
        
      }
      
      ITIS_id<-tsn_n$tsn
      ITIS_name<-tsn_n$scientificName
      ITIS_is_valid<-tsn_n$nameUsage
      
      # Get the upstream taxonomic information
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      j<-NULL
      t_13<-1
      
      while( is.null(j) && t_5 <= n_times){
        
        try(j<-itis_hierarchy(tsn_n$tsn,what="full",silent=TRUE))
        t_13 <- t_13 + 1
      }
      #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
      
      ITIS_Phylum<-toupper(as.character(j[j$rankname=="phylum",4])) # Change from lower to upper chase
      ITIS_Class<-toupper(as.character(j[j$rankname=="class",4]))
      ITIS_Order<-toupper(as.character(j[j$rankname=="order",4]))
      ITIS_Family<-toupper(as.character(j[j$rankname=="family",4]))
      
    }
  }
  
  # b.3.2 Unify the taxonomic information from ITIS---- 
  ITIS_data<-data.frame(ITIS_Present,ITIS_is_valid,ITIS_id,ITIS_name,
                        ITIS_Phylum,ITIS_N_syn,ITIS_syn,
                        ITIS_Class,ITIS_Order,ITIS_Family)
  
  # Should we retrieve synonim information from GBIF?
  if(Gbif==TRUE){
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    key_1 <- NULL
    t_6 <- 1
    
    while(is.null(key_1) && t_6 <= n_times){
      
      try(key_1 <- get_gbifid_(sci=spp.x)[[1]],silent = TRUE) # get the taxon key
      t_6 <- t_6 + 1
    }
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    
    if (length(key_1)==0){
      
      GBIF_Present<-"No"
      GBIF_id<-NA
      
      GBIF_name<-NA
      
      GBIF_Phylum<-NA
      GBIF_Class<-NA
      GBIF_Order<-NA
      GBIF_Family<-NA
      
      GBIF_syn<-NA
      GBIF_N_syn<-NA
      
    } else{
      
      GBIF_id<-paste(key_1[key_1$status=="ACCEPTED" & key_1$matchtype=="EXACT",20],collapse="-")
      GBIF_name<-paste(unique(key_1[key_1$specieskey==GBIF_id,13]),collapse="-")
      
      GBIF_Phylum<-toupper(paste(unique(key_1[key_1$specieskey==GBIF_id,9]),collapse="-"))
      GBIF_Class<-toupper(paste(unique(key_1[key_1$specieskey==GBIF_id,22]),collapse="-"))
      GBIF_Order<-toupper(paste(unique(key_1[key_1$specieskey==GBIF_id,10]),collapse="-"))
      GBIF_Family<-toupper(paste(unique(key_1[key_1$specieskey==GBIF_id,11]),collapse="-"))
      
      if (is.na(GBIF_name)){
        
        GBIF_Present<-"No"
        GBIF_id<-NA
        
        GBIF_name<-NA
        
        GBIF_Phylum<-NA
        GBIF_Class<-NA
        GBIF_Order<-NA
        GBIF_Family<-NA
        
        GBIF_syn<-NA
        GBIF_N_syn<-NA
        
      } 
      
      if (GBIF_name==""){
        
        GBIF_Present<-"No"
        GBIF_id<-NA
        
        GBIF_name<-NA
        
        GBIF_Phylum<-NA
        GBIF_Class<-NA
        GBIF_Order<-NA
        GBIF_Family<-NA
        
        GBIF_syn<-NA
        GBIF_N_syn<-NA
        
      } else {
        
        GBIF_Present<-"Yes" 
        GBIF_syn<-paste(key_1[key_1$status=="SYNONYM",13],collapse=";")
        if(length(key_1[key_1$status=="SYNONYM",1])==0){
          GBIF_N_syn<-NA
        }else{
          GBIF_N_syn<-length(key_1[key_1$status=="SYNONYM",1]) 
        }
      }
    }
    
    GBif_data<-data.frame(GBIF_Present,GBIF_id,GBIF_name,
                          GBIF_N_syn,GBIF_syn,GBIF_Phylum,
                          GBIF_Class,GBIF_Order,GBIF_Family)
  }
  
  # C. Return the Taxonomic information
  if(Gbif==TRUE){
    Tax_dat<-cbind(Or_name=spp.x,IUCN_data,ITIS_data,GBif_data)
    Spp_syn<-c(spp.x,Tax_dat[,colnames(Tax_dat) %in% c("IUCN_name","IUCN_syn","ITIS_name",
                                                       "ITIS_syn","GBIF_name","GBIF_syn")]) %>% paste(collapse = ";") %>% 
      strsplit(split = ";") %>% unlist()
    
    Spp_syn<-Spp_syn[-c(Spp_syn%>%grep(pattern = "NA"))]
    Spp_syn<-Spp_syn[!duplicated(Spp_syn)]
    
    return(list(Spp_syn=Spp_syn,
                IUCN_spp=IUCN_data$IUCN_name,
                Or_name=spp.x,
                TaxDat=Tax_dat))
  }else{
    Tax_dat<-cbind(Or_name=spp.x,IUCN_data,ITIS_data)
    Spp_syn<-c(spp.x,Tax_dat[,colnames(Tax_dat) %in% c("IUCN_name","IUCN_syn","ITIS_name",
                                                       "ITIS_syn","GBIF_name","GBIF_syn")]) %>% paste(collapse = ";") %>% 
      strsplit(split = ";") %>% unlist()
    
    Spp_syn<-Spp_syn[-c(Spp_syn%>%grep(pattern = "NA"))]
    Spp_syn<-Spp_syn[!duplicated(Spp_syn)]
    
    return(list(Spp_syn=Spp_syn,
                IUCN_spp=IUCN_data$IUCN_name,
                Or_name=spp.x,
                TaxDat=Tax_dat))
    
  }
  # Message 
  print("Taxonomic search done")
}

# End of the script