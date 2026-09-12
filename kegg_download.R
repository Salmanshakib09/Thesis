library(msigdbr)
library(readr)
kegg= msigdbr(db_species = "HS", species = "Homo sapiens", 
              collection ="C2", subcollection ="CP:KEGG_MEDICUS")
unique(kegg$db_version)
dir.create("data/pathways",
           recursive = TRUE, 
           showWarnings = FALSE)
write_csv(kegg,
          "data/pathways/kegg_medicus.csv")