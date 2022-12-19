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

palette_km <- function(group_names){
  
  groups <- length(unique(group_names))
  cols<-c("green", "#009966", "Blue", "Purple", "Red")
  
  
  if(groups==1){
    return(cols[1])
  }
  else if(groups==2){
    return(cols[c(1,5)])
  }
  else if(groups==3){
    return(cols[c(1,3,5)])
  }
  else if(groups==4){
    return(cols[c(1,2,4,5)])
  }
  else if(groups==5){
    return(cols)
  }
}

KM_curve <- function(mRNA_counts,clinical_data, gene="COL7A1", gene_notation="SYMBOL", 
                     grouping="equal range", time_unit="years",histogram_bins=20,
                     group_names=c("1:Low", "2:Low to Medium", "3:Medium",
                                   "4:Medium to High", "5:High"), 
                     censoring=T, facet_by=FALSE, use_outliers=T,
                     ranges, histogram_limits_x=NULL,
                     histogram_limits_y=NULL, survival_limit_x=NULL, only_km = F){
  
  require(tidyverse)
  require(survival)
  require(survminer)
  require(ggsci)
  ranges_hist=NULL
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
      
      ranges_hist <- levels(cut(as_vector(KM_mRNA[,gene]), 5*group_number))
      ranges_hist <- sort(unique(as.double(str_replace(str_split_fixed(ranges_hist, "\\,",2 ), "[\\(\\]]", ""))))
      ranges_hist <- ranges_hist[3:length(ranges_hist)-1]
    }
    else if (use_outliers==F){
      KM_no_outliers <- KM_mRNA %>% 
        select(gene) %>% 
        filter(!(abs(get(gene) - median(get(gene))) > 6*mad(get(gene))))
      
      ranges <- levels(cut(as_vector(KM_no_outliers[,gene]), group_number)) 
      ranges <- as.double(str_sub(str_split_fixed(ranges, "\\,",2 )[-length(ranges),2], 1, -2))
      ranges <- c(-Inf, ranges, +Inf)
      
      
      ranges_hist <- levels(cut(as_vector(KM_no_outliers[,gene]), 5*group_number))
      ranges_hist <- sort(unique(as.double(str_replace(str_split_fixed(ranges_hist, "\\,",2 ), "[\\(\\]]", ""))))
      ranges_hist <- ranges_hist[3:length(ranges_hist)-1]
      
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
    
    KM_no_outliers <- KM_mRNA %>% 
      select(gene) %>% 
      filter(!(abs(get(gene) - median(get(gene))) > 6*mad(get(gene))))
    
    ranges<-c(-Inf, ranges, Inf)
    
    KM_mRNA <- KM_mRNA %>% 
      mutate("expression"= cut(x = as_vector(KM_mRNA[,gene]), 
                               breaks = ranges,labels=group_names)) %>% 
      mutate("expression" = as.factor(expression))
  }
  
  else{stop('Value of "grouping"  must be either "equal range","set range" or "equal groups".')}
  
  breakpoints_hist=NULL
  if(length(ranges)==3 & grouping=="set range"){
    range_expression <- range(as.vector(KM_mRNA[,gene]))
    breakpoints_hist<-seq(range_expression[1],(range_expression[2]+ranges[2]), ((ranges[2]-range_expression[1])/histogram_bins))
    
    KM_hist <- ggplot(KM_mRNA, aes(x=get(gene), fill=expression))+
      geom_histogram(position='identity', 
                     breaks=breakpoints_hist,
                     alpha=0.7, color="Black")+
      theme_minimal()+
      scale_fill_manual(values = palette_km(group_names))+
      xlab(gene)+
      coord_cartesian(xlim = histogram_limits_x,ylim = histogram_limits_y, default = FALSE, expand = T)
  }
  else{
    KM_hist <- ggplot(KM_mRNA, aes(x=get(gene), fill=expression))+
      geom_histogram(position='identity', 
                     breaks=ranges_hist,
                     alpha=0.7, color="Black")+
      theme_minimal()+
      scale_fill_manual(values = palette_km(group_names))+
      xlab(gene)+
      coord_cartesian(xlim = histogram_limits_x,ylim = histogram_limits_y, default = FALSE, expand = T)
  }
  
  
  if (facet_by==FALSE){
    KM_fit <- survfit(Surv(survival_time, vital_status) ~ expression, data = KM_mRNA)
    
    KM_plot <- ggsurvplot(KM_fit, xlab = str_to_sentence(time_unit), 
                          ylab = "Overall survival probability", pval = T,risk.table = T,
                          legend="none", censor=censoring, tables.height=0.35, data=KM_mRNA,
                          palette=palette_km(group_names), xlim=survival_limit_x)
    if(!only_km){
      return(list(KM_hist, KM_plot))
    }
    else{return(KM_plot)
    }
  }
  else{
    
    KM_fit <- survfit(Surv(survival_time, vital_status) ~ expression, data = KM_mRNA)
    
    KM_plot <- ggsurvplot(fit = KM_fit, xlab =  str_to_sentence(time_unit), 
                          ylab = "Overall survival probability", 
                          censor=censoring,  data=KM_mRNA, facet.by = facet_by,
                          palette=palette_km(group_names), add.all = TRUE)
    
    form<-formula(paste0("Surv(survival_time, vital_status) ~ expression +", facet_by))
    
    KM_table <- survfit(form, data = KM_mRNA)
    
    KM_table<-enframe(KM_table[[10]]) %>% 
      mutate("expression"=str_split_fixed(name,"\\,",n=2)[,1]) %>% 
      mutate("expression"=as.factor(str_split_fixed(expression,"\\=",n=2)[,2])) %>% 
      mutate("strata"=str_split_fixed(name,"\\,",n=2)[,2]) %>% 
      mutate("strata"=str_split_fixed(strata,"\\=",n=2)[,2]) %>% 
      arrange(expression) %>% 
      select(-name) %>% 
      ggplot()+
      geom_text(aes(y=reorder(expression, desc(expression)), x=strata, label=value), size=6)+
      theme_bw() + 
      theme(panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            strip.background = element_rect(fill="white", size = 1),
            axis.text.x = element_text(color="black", size = 12),
            axis.text.y = element_text(color="black", size = 12),
            axis.title = element_text(color="black", size = 14),
            strip.text = element_text(color="black", size = 12))+
      ggtitle("Numbers of patients") +
      xlab(str_to_title(facet_by)) + ylab("Expression Levels")+
      scale_x_discrete(position = "top")+
      scale_y_discrete(limit=levels(expression))
    
    KM_mRNA<-KM_mRNA %>% filter(!is.na(get(facet_by)))
    
    l<-count(KM_mRNA, holder=get(facet_by))
    colnames(l)[1]=facet_by
    
    
    KM_hist <- ggplot(KM_mRNA, aes(x=get(gene), fill=expression))+
      geom_histogram(position='identity',
                     breaks=breakpoints_hist,
                     alpha=0.7, 
                     mapping = aes(y=..density..),
                     color="Black")+
      theme_minimal()+
      scale_fill_manual(values = palette_km(group_names))+
      xlab(gene)+
      coord_cartesian(xlim = histogram_limits_x,ylim = histogram_limits_y, default = FALSE, expand = T)+
      facet_wrap(~get(facet_by), nrow = 1)+
      geom_text(data= l, 
                aes(x=(mean(histogram_limits_x)), y=histogram_limits_y[2],
                    label=paste("n= ",n)), size=5,nudge_x = 0, inherit.aes=FALSE )
    
    return(list(KM_hist,KM_plot,KM_table))
    
  }
  
}

