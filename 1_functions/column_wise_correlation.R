## description:
# Function can take two datasets and compute correlation coefficient between all columns.
# Alternatively, one dataset can be the input, and genes will be correlated with each other.
# There is also an option to just correlate one specific gene with rest of genes in dataset. 
# For later two cases, "set2" argument should be declared as "self" or "gene" resp.
# Data can be also be normalized by function, and prior count can be changed if log2 normalization is selected. 
# Method for calculation of correlation can be changed. 
# Function can additionally return histogram of correlation coefficients, with overlaying normal curve.
# Preferably, input is dataset generated from "read data" functions from this folder. 
#
# note: normalize_cpm_log function is required to perform normalization inside a function.
##................................................................................................................

column_wise_correlation <- function(set1, set2="gene", gene="COL7A1", normalization="none", 
                                    prior.count=0.000001, method="pearson", return="data", binwidth_hist=0.02){
  require(tidyverse)
  require(ggpubr)
  require(psych) 

  if(typeof(set2) == "character"){
    if(set2 == "self"){
      self_remove=TRUE
      
      if(normalization == "log"){
        set1_norm <- normalize_cpm_log(set1,log.2=TRUE,prior.count=prior.count)
      }
      else if (normalization == "cpm"){
        set1_norm <- normalize_cpm_log(set1)
      }
      else if ((normalization == "none")){
        set1_norm <- set1
      }
      else{stop('Mode of normalization must be "log" or "cpm". If data is already normalized, use value "none".')
      }
      set2_norm <- set1_norm
    }
    else if (set2 == "gene"){
      self_remove=TRUE
      if(normalization == "log"){
        set1_norm <- normalize_cpm_log(set1,log.2=TRUE,prior.count=prior.count)
      }
      else if (normalization == "cpm"){
        set1_norm <- normalize_cpm_log(set1)
      }
      else if ((normalization == "none")){
        set1_norm <- set1
      }
      else{stop('Mode of normalization must be "log" or "cpm". If data is already normalized, use value "none".')
      }
      
      named <- names(set1_norm)
      names(named) <- str_split_fixed(names(set1_norm),"\\|",2)[,1]
      set2_norm <- set1_norm %>% 
        select(1, named[gene])
    }
  }
  else if(typeof(set2) == "list"){
    self_remove=FALSE
    if(normalization == "log"){
      set1_norm <- normalize_cpm_log(set1,log.2=TRUE,prior.count=prior.count)
      set2_norm <- normalize_cpm_log(set2,log.2=TRUE,prior.count=prior.count)
      
    }
    else if (normalization == "cpm"){
      set1_norm <- normalize_cpm_log(set1)
      set2_norm <- normalize_cpm_log(set2)
    }
    else if (normalization == "none"){
      set1_norm <- set1
      set2_norm <- set2
    }
    else{stop('Mode of normalization must be "log" or "cpm". If data is already normalized, use value "none".')
    }
    set1_norm <- set1_norm %>% 
      semi_join(set2_norm, by="patient")
    set2_norm <- set2_norm %>% 
      semi_join(set1_norm, by="patient")
  }
  
  correlated_sets <- corr.test(x =set1_norm[,-1] ,y=set2_norm[,-1],  adjust="BH", alpha=0.05, method = method)
  
  correlated_sets_df<-as_tibble(correlated_sets$r)
  
  correlated_sets_df <- correlated_sets_df %>% 
    mutate("set1_genes"=rownames(correlated_sets$r)) %>% 
    select(set1_genes, everything()) %>% 
    pivot_longer(cols= 2:(ncol(correlated_sets_df)+1), names_to="set2_genes",
                 values_to="correlation") %>% 
    filter(!is.na(correlation))
  
  correlated_sets_p<-as_tibble(correlated_sets$p)
  
  correlated_sets_p <- correlated_sets_p %>% 
    mutate("set1_genes"=rownames(correlated_sets$p)) %>% 
    select(set1_genes, everything()) %>% 
    pivot_longer(cols= 2:(ncol(correlated_sets_p)+1), names_to="set2_genes",
                 values_to="p.adj")
  
  correlated_sets_df <- correlated_sets_df %>%
    left_join(correlated_sets_p, by=c("set1_genes", "set2_genes"))
  
  
  if(self_remove == TRUE){
    correlated_sets_df <- correlated_sets_df %>% 
      mutate("self_rem"=str_split_fixed(set1_genes,"\\|",2)[,1])%>% 
      filter(!set1_genes == set2_genes) %>% 
      filter(!self_rem == set2_genes) %>% 
      select(-self_rem)
    
  }
  if(return == "data"){
    return(correlated_sets_df)
  }
  else if (return == "histogram")
    ploted_histogram <-correlated_sets_df %>% 
    group_by(set2_genes) %>% 
    nest(data =c(set1_genes, p.adj, correlation)) %>% 
    mutate(y=map(data, ~ dnorm(
      .$correlation, mean=median(.$correlation), sd=mad(.$correlation)
    ) * binwidth_hist * sum(!is.na(.$correlation)))) %>% 
    unnest(c(data,y)) 
  
  ploted_histogram <- ggplot(data=ploted_histogram, aes(x=correlation))+
    geom_histogram(binwidth=binwidth_hist, fill="#33CCFF", color="#000000")+
    geom_line(aes(y=y), color="#FF0000", lwd=1.5, alpha=0.6)+
    labs(title=paste0("Gaussian distribution: mu=",
                      signif(mean(ploted_histogram$correlation), 3),
                      ", sigma=", signif(mad(ploted_histogram$correlation), 3)))+
    xlab("Correlation coefficient")+
    theme_pubr()
  return(ploted_histogram)
}
