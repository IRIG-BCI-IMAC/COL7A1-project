# 95% confidence interval formula
CI_median_95 <- function(gene){
low <- round(length(gene)*0.5 - 1.96*sqrt(length(gene)*0.5*0.5))
high <- round(length(gene)*0.5 + 1.96*sqrt(length(gene)*0.5*0.5))

low <- sort(gene)[low]
high <- sort(gene)[high]

return(c(low, high))
}

CI_median_95(subset(hist_data, Gene=="ENSG00000114270.17")$Expression)
CI_median_95(subset(hist_data, Gene!="ENSG00000114270.17")$Expression)
median(subset(hist_data, Gene!="ENSG00000114270.17")$Expression)
