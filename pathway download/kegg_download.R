library(tidyverse)
library(msigdbr)
kegg= msigdbr(db_species = "HS", 
              species = "Homo sapiens", 
              collection = "C2",
              subcollection = "CP:KEGG_MEDICUS"
              )
kegg_gmt= kegg |> 
          select(gs_name, gs_description, gene_symbol)|>
          distinct()|>
          group_by(gs_name, gs_description)|>
          summarise(genes= paste(gene_symbol, collapse = "\t"),
                    .groups = "drop"
          )|>
          mutate(
            line=paste(gs_name, gs_description, genes, sep = "/t")
          )
dir.create("data/pathways",
           recursive=TRUE,
           showWarnings = FALSE)
write_lines( kegg_gmt$line,
            "data/pathways/kegg_medicus.gmt")