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

```R
listOfGwasIds <- read.table("TRAITS.txt", header = T)
for(i in 1:nrow(listOfGwasIds))
{
  instrumentId <- as.character(listOfGwasIds$id[i])
  instrumentName <- as.character(listOfGwasIds$study[i])
  tag <- paste("INSTRUMENT IS ",instrumentId," AT i = ",i, sep = "")
  print(tag)
  flagged <- "nope"
  Exp_data <- extract_instruments(outcomes=instrumentId, p1 = 5e-08, clump = TRUE, p2 = 5e-08,
                                  r2 = 0.001, kb = 10000, access_token = token,
                                  force_server = TRUE)
  skip <- ifelse(length(Exp_data$beta.exposure) < 10, 1, 0)
  if(skip == 0)
  {
    ncase.exposure <- listOfGwasIds$ncase[i]
    ncontrol.exposure <- listOfGwasIds$ncontrol[i]
    dat <- harmonise_data(exposure_dat=Exp_data, outcome_dat=oc_meta, action=2)
    dat$samplesize.exposure <- listOfGwasIds$samplesize[i]
    if(ncase.exposure > 0) 
    {
      dat$ncase.exposure <- ncase.exposure
      dat$ncontrol.exposure <- ncontrol.exposure
      dat1<-subset(dat, dat$eaf.exposure!="NA")
      dat1$r.exposure<- get_r_from_lor(dat1$beta.exposure, dat1$eaf.exposure, dat1$ncase.exposure, dat1$ncontrol.exposure, 0.01,  model = "logit")
      } else {
      dat1<- dat
      dat1$r.exposure<- get_r_from_pn(dat1$pval.exposure, dat1$samplesize.exposure)
    }
    steiger <- steiger_filtering(dat1)
    sig<-subset(steiger, steiger$steiger_dir==TRUE)
    presso <-mr_presso(BetaOutcome = "beta.outcome", BetaExposure = "beta.exposure", SdOutcome = "se.outcome", SdExposure = "se.exposure", OUTLIERtest = TRUE, DISTORTIONtest = TRUE, data = sig, NbDistribution = 1000,  SignifThreshold = 0.05)
    capture.output(print(presso), file = paste(instrumentName, "presso.txt", sep = ""))
    res <- mr(sig)
    print(res)
    res_single <- mr_singlesnp(sig)
    het <- mr_heterogeneity(dat)
    print(het)
    ple <- mr_pleiotropy_test(dat)
    print(ple)
    write.table(res, file = paste(instrumentName,"_1res.txt",sep = ""), quote = F, sep = ",")
    write.table(het, file = paste(instrumentName,"_2het.txt",sep = ""), quote = F, sep = ",")
    write.table(ple, file = paste(instrumentName,"_3ple.txt",sep = ""), quote = F, sep = ",")
    res_single <- mr_singlesnp(sig)
    p5 <- mr_forest_plot(res_single)
    ggsave(p5[[1]], file= paste(instrumentName,"_forest.jpg",sep=""), width=7, height=12)
    write.table(res_single, file = paste(instrumentName,"_4res_single.txt",sep = ""), quote = F, sep = ",")
    out <- directionality_test(sig)
    write.table(out, file = paste(instrumentName,"_5dir.txt",sep = ""), quote = F, sep = ",")
  }
  else
  {
    print("FAIL")
  }
}
```
