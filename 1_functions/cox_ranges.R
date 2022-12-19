## description: 
# Function used for grouping patients according to gene expression levels:
#   grouping: could be done in 3 ways: 
#     "equal range" cuts total range of expression in up to 10 different levels
#     "equal groups" separates total number of patients in up to 10 equal groups
#     "set range" allows for the user to appoint desired cut-off points for group division.
#
#   use_outliers: allows user to decide where outlier's would be used in calculation
# of total range of expression. This is usefull when "equal range" is grouping option.
#
#   facet_by: allows user to facet Kaplan-Meier curve according to certain characteristics
# found in clinical data (example: stage).
#
#   Current pallete is deffined in function "palette_km", and can be changed.
##..............................................


cox_ranges <- function(mRNA_counts,clinical_data, gene="COL7A1", gene_notation="SYMBOL", 
                     grouping="equal range",
                     group_names=c("1:Low", "2:Low to Medium", "3:Medium",
                                   "4:Medium to High", "5:High"),  use_outliers=T){
  
  require(tidyverse)
  require(survival)
  require(survminer)
  require(ggsci)
  
  if(gene_notation=="SYMBOL"){
    names(mRNA_counts) <- str_split_fixed(names(mRNA_counts),"\\|",2)[,1]
  }
  else if(gene_notation=="ENTREZ"){
    names(mRNA_counts) <- c("patient", str_split_fixed(names(mRNA_counts),"\\|",2)[-1,2])
  }
  else{stop('Currently suported gene notations are"SYMBOL" and "ENTREZ"')}
  
  KM_mRNA <- mRNA_counts %>% 
    select(patient, sym(gene)) %>% 
    inner_join(clinical_data, by="patient")


  group_number <- length(group_names)
  
  if(grouping=="equal range"){
    if (use_outliers==T){
    KM_mRNA <- KM_mRNA %>% 
    mutate("expression"= cut(x = as_vector(KM_mRNA[,gene]), 
                             breaks = group_number, labels=group_names)) %>% 
      mutate("expression" = as.factor(expression))
    }
    else if (use_outliers==F){
      KM_no_outliers <- KM_mRNA %>% 
        select(gene) %>% 
        filter(!(abs(get(gene) - median(get(gene))) > 6*mad(get(gene))))
      
      ranges <- levels(cut(as_vector(KM_no_outliers[,gene]), group_number)) 
      ranges <- as.double(str_sub(str_split_fixed(ranges, "\\,",2 )[-length(ranges),2], 1, -2))
      ranges <- c(-Inf, ranges, +Inf)
      
      KM_mRNA <- KM_mRNA %>% 
        mutate("expression"= cut(x = as_vector(KM_mRNA[,gene]), 
                                 breaks = ranges, labels=group_names)) %>% 
        mutate("expression" = as.factor(expression))
    }
  }
  else if(grouping=="equal groups"){
    KM_mRNA <- KM_mRNA %>% 
      arrange(get(gene)) %>% 
      mutate("expression"=ceiling((row_number()*group_number)/nrow(KM_mRNA))) %>% 
      mutate("expression"=case_when(expression==1~group_names[1], expression==2~group_names[2], 
                                  expression==3 ~group_names[3],expression==4 ~group_names[4],
                                  expression==5 ~group_names[5],expression==6 ~group_names[6],
                                  expression==7 ~group_names[7],expression==8 ~group_names[8],
                                  expression==9 ~group_names[9],expression==10 ~group_names[10])) %>% 
      mutate("expression" = as.factor(expression))
  }
  
  else if(grouping=="set range"){
    
    if(length(group_names)!=(length(ranges)+1)){stop("for argument 'ranges' provide one less value than for 'group names'")}
    
    ranges<-c(-Inf, ranges, Inf)
    
    KM_mRNA <- KM_mRNA %>% 
      mutate("expression"= cut(x = as_vector(KM_mRNA[,gene]), 
                               breaks = ranges,labels=group_names)) %>% 
      mutate("expression" = as.factor(expression))
  }
  
  else{stop('Value of "grouping"  must be either "equal range","set range" or "equal groups".')}
  
    cox_res <- coxph(Surv(time = survival_time, event = vital_status)~expression, data = KM_mRNA[,-1])
    
    return(cox_res)
}
    
    
    
    
    