# RBD_MR

## Set up. 
```R
devtools::install_github("MRCIEU/TwoSampleMR")
devtools::install_github("MRCIEU/MRInstruments")

require(TwoSampleMR)
require(ggplot2)
require(devtools)
require(MRInstruments) 
require(data.table)

token <- ieugwasr::check_access_token() # MRBase requires a token to use. 
````
## Preparing outcomes
**If you need to check available summary stats in MRBase, I recommend doing the following:**
```R
ao = available_outcomes()
write.table(ao, file="MR_available_outcomes.txt", col.names=T, row.names=F, sep="\t", quote=F)
````

From MR_available_outcomes.txt, find the traits you are testing and their corresponding ID. Save these IDs to TRAITS.txt.  
*(Many of these summary stats have the same name and PubMed ID as the LDHub summary stats. This makes it easier if you are testing genetically correlated traits, as I did.)*   

## Run MR
This loop is an adapation from Sara Bandres-Ciga (www.github.com/sarabandres).  
*Right now it still has bugs. Working on it.* 

```R
listOfGwasIds <- read.table("TRAITS.txt", header = T)
for(i in 1:nrow(listOfGwasIds))
{
  instrumentId <- as.character(listOfGwasIds$id[i])
  tag <- paste("INSTRUMENT IS ",instrumentId," AT i = ",i, sep = "")
  print(tag)
  flagged <- "nope"
  Exp_data <- extract_instruments(outcomes=instrumentId, p1 = 5e-08, clump = TRUE, p2 = 5e-08,
                                  r2 = 0.001, kb = 10000, access_token = token,
                                  force_server = TRUE)
  skip <- ifelse(length(Exp_data$beta.exposure) < 1, 1, 0)
  if(skip == 0)
  {
    dat <- harmonise_data(exposure_dat=Exp_data, outcome_dat=Out_data, action=2)
    res <-mr(dat)
    print(res)
    het <- mr_heterogeneity(dat)
    print(het)
    ple <- mr_pleiotropy_test(dat)
    print(ple)
    write.table(res, file = paste(instrumentId,"_1res.txt",sep = ""), quote = F, sep = ",")
    write.table(het, file = paste(instrumentId,"_2het.txt",sep = ""), quote = F, sep = ",")
    write.table(ple, file = paste(instrumentId,"_3ple.txt",sep = ""), quote = F, sep = ",")
    res_single <- mr_singlesnp(res)
    p1 <- mr_forest_plot(res_single)
    ggsave(p1[[1]], file=paste(instrumentId,"_forest.png", width=7, height=7)
    write.table(res_single, file = paste(instrumentId,"_4res_single.txt",sep = ""), quote = F, sep = ",")
    out <- directionality_test(dat)
  }
  else
  {
    print("FAIL")
  }
}
```
