library(tidyverse)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(clusterProfiler)


# Import data -------------------------------------------------------------

files <-
  list.files(path = "Data/Spots_Matheus_L_nov24/",
           pattern = "protein.csv",full.names = T,recursive = T)


tbl <- 
  map(files, read_csv2) %>% 
  map(janitor::clean_names)


membranas_ves <-
  tbl[9:10] %>% 
  bind_rows()

citosol_ves <-
  tbl[11:12] %>% 
  bind_rows()

# Proteinas ---------------------------------------------------------------

memprot <- membranas_ves$protein_accession
citosolprot <- citosol_ves$protein_accession

AnnotationDbi::columns(org.Hs.eg.db)

memprot_tbl <- select(org.Hs.eg.db,
                      keys = memprot,
                      columns = c("SYMBOL","ENTREZID","UNIPROT"),
                      keytype = "UNIPROT"
                      )
citosolprot_tbl <- AnnotationDbi::select(org.Hs.eg.db,
                      keys = citosolprot,
                      columns = c("SYMBOL","ENTREZID","UNIPROT"),
                      keytype = "UNIPROT")

# GO ----------------------------------------------------------------------

memprot_go <- enrichGO(gene = unique(memprot_tbl$SYMBOL),
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "ALL",
                       minGSSize = 1,
                       maxGSSize = 500)

citosolprot_go <- enrichGO(gene = unique(citosolprot_tbl$SYMBOL),
                       OrgDb = org.Hs.eg.db,
                       keyType = "SYMBOL",
                       ont = "ALL",
                       minGSSize = 1,
                       maxGSSize = 500)


# STRING ------------------------------------------------------------------

clipr::write_clip(unique(memprot_tbl$SYMBOL))
clipr::write_clip(unique(citosolprot_tbl$SYMBOL))
