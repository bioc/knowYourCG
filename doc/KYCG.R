## ----ky1, load-depenencies, results="hide", message=FALSE, warning=FALSE------
library(knowYourCG)
library(sesameData)
sesameDataCache(data_titles=c("genomeInfo.hg38","genomeInfo.mm10",
                  "KYCG.MM285.tissueSignature.20211211",
                  "MM285.tissueSignature","MM285.address",
                  "probeIDSignature","KYCG.MM285.designGroup.20210210",
                  "KYCG.MM285.chromHMM.20210210",
                  "KYCG.MM285.TFBSconsensus.20220116",
                  "KYCG.MM285.HMconsensus.20220116",
                  "KYCG.MM285.chromosome.mm10.20210630"
                  ))

## ----ky2, message=FALSE-------------------------------------------------------
query <- getDBs("MM285.designGroup")[["PGCMeth"]]
head(query)

## ----ky3, fig.width=8, fig.height=5, message=FALSE----------------------------
dbs <- c("KYCG.MM285.chromHMM.20210210",
         "KYCG.HM450.TFBSconsensus.20211013",
         "KYCG.MM285.HMconsensus.20220116",
         "KYCG.MM285.tissueSignature.20211211",
         "KYCG.MM285.chromosome.mm10.20210630",
         "KYCG.MM285.designGroup.20210210")
results_pgc <- testEnrichment(query,databases = dbs,platform="MM285")
head(results_pgc)

## ----ky9, echo = FALSE, results="asis"----------------------------------------
library(knitr)
df = data.frame(
    c("Correlation-based","Set Enrichment Analysis"),
    c("Set Enrichment Analysis","Fisher's Exact Test")
)
colnames(df) <- c("Continuous Database Set", "Discrete Database Set")
rownames(df) <- c("Continuous Query", "Discrete Query")
kable(df, caption="Four knowYourCG Testing Scenarios")

## ----ky10, run-test-single, echo=TRUE, eval=TRUE, message=FALSE---------------
library(SummarizedExperiment)

## prepare a query
df <- rowData(sesameDataGet('MM285.tissueSignature'))
query <- df$Probe_ID[df$branch == "fetal_brain" & df$type == "Hypo"]

results <- testEnrichment(query, "TFBS", platform="MM285")
results %>% dplyr::filter(overlap>10) %>% head

## prepare another query
query <- df$Probe_ID[df$branch == "fetal_liver" & df$type == "Hypo"]
results <- testEnrichment(query, "TFBS", platform="MM285")
results %>% dplyr::filter(overlap>10) %>%
    dplyr::select(dbname, estimate, test, FDR) %>% head

## ----ky5, list-data, eval=TRUE, echo=TRUE-------------------------------------
listDBGroups("MM285")

## ----ky6, cache-data, eval=TRUE, warning=FALSE--------------------------------
dbs <- getDBs("MM285.design")

## ----ky7, view-data1, eval=TRUE, warning=FALSE--------------------------------
str(dbs[["PGCMeth"]])

## ----ky8, message=FALSE-------------------------------------------------------
df <- rowData(sesameDataGet('MM285.tissueSignature'))
query <- df$Probe_ID[df$branch == "B_cell"]
head(query)

## ----ky16, fig.width=7, fig.height=6, echo=TRUE, warning=FALSE, message=FALSE----
query <- names(sesameData_getProbesByGene("Dnmt3a", "MM285"))
results <- testEnrichment(query, 
    buildGeneDBs(query, max_distance=100000, platform="MM285"),
    platform="MM285")
main_stats <- c("dbname","estimate","gene_name","FDR", "nQ", "nD", "overlap")
results[,main_stats]

## ----ky28, warning=FALSE, message=FALSE, fig.width=5, fig.height=4------------
query <- names(sesameData_getProbesByGene("Dnmt3a", "MM285"))
dbs <- c("KYCG.MM285.chromHMM.20210210","KYCG.HM450.TFBSconsensus.20211013",
         "KYCG.MM285.chromosome.mm10.20210630")
results <- testEnrichment(query,databases=dbs,
                          platform="MM285",include_genes=TRUE)
main_stats <- c("dbname","estimate","gene_name","FDR", "nQ", "nD", "overlap")
results[,main_stats] %>% 
    head()

## ----ky18, message=FALSE------------------------------------------------------
library(gprofiler2)
df <- rowData(sesameDataGet('MM285.tissueSignature'))
query <- df$Probe_ID[df$branch == "fetal_liver" & df$type == "Hypo"]
res <- testGO(query, platform="MM285",organism = "mmusculus")
head(res$result)

## ----ky29, eval = TRUE,message=FALSE------------------------------------------
df <- rowData(sesameDataGet('MM285.tissueSignature'))
probes <- df$Probe_ID[df$branch == "fetal_liver" & df$type == "Hypo"]
res <- testProbeProximity(probeIDs=probes)
head(res)

## ----ky21, run-test-data, echo=TRUE, eval=TRUE, message=FALSE-----------------
query <- getDBs("KYCG.MM285.designGroup")[["TSS"]]

## ----ky22, echo=TRUE, eval=TRUE, message=FALSE--------------------------------
sesameDataCache(data_titles = c("KYCG.MM285.seqContextN.20210630"))
res <- testEnrichmentSEA(query, "MM285.seqContextN")
main_stats <- c("dbname", "test", "estimate", "FDR", "nQ", "nD", "overlap")
res[,main_stats]

## ----ky23, warning=FALSE, eval=TRUE,message=FALSE-----------------------------
library(sesame)
sesameDataCache(data_titles = c("MM285.1.SigDF"))
beta_values <- getBetas(sesameDataGet("MM285.1.SigDF"))
res <- testEnrichmentSEA(beta_values, "MM285.chromHMM")
main_stats <- c("dbname", "test", "estimate", "FDR", "nQ", "nD", "overlap")
res[,main_stats] 

## -----------------------------------------------------------------------------
sessionInfo()

