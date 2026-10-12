source("analysis/00_setup.R")

# Standard GMT: pathway name, description, then gene symbols, separated by tabs.
write_reactome_gmt <- function(db, path="data/pathways/reactome.gmt") {
  rows <- split(seq_len(nrow(db)), db$gs_name)
  lines <- vapply(names(rows), function(name) {
    i <- rows[[name]]
    description <- if ("gs_description" %in% names(db)) db$gs_description[i[1]] else NA_character_
    if (is.na(description) || !nzchar(description)) {
      description <- if ("gs_url" %in% names(db)) db$gs_url[i[1]] else "na"
    }
    if (is.na(description) || !nzchar(description)) description <- "na"
    description <- gsub("[\t\r\n]", " ", description)
    genes <- unique(db$gene_symbol[i])
    if (anyNA(genes) || any(!nzchar(genes)) || any(grepl("[\t\r\n]", genes))) stop("Invalid GMT gene symbols")
    paste(c(name, description, genes), collapse="\t")
  }, character(1))
  dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
  writeLines(lines, path, useBytes=TRUE)
  invisible(path)
}

load_reactome <- function() {
  path <- 'data/raw/reactome-msigdbr.rds'
  if (!file.exists(path)) {
    db <- msigdbr::msigdbr(species='Homo sapiens',collection='C2',subcollection='CP:REACTOME')
    dir.create(dirname(path), recursive=TRUE, showWarnings=FALSE)
    saveRDS(db,path)
  } else db <- readRDS(path)
  write_reactome_gmt(db)
  list(sets=lapply(split(db$gene_symbol,db$gs_name),unique),
       version=unique(db$db_version))
}

reactome_stage <- function(directory=stage_dir) {
  input <- read_stage("input", directory)
  db <- load_reactome()
  sets <- lapply(db$sets, intersect, rownames(input$expression))
  eligible <- sort(names(sets)[lengths(sets)>=20 & lengths(sets)<=200])
  if (length(eligible)<10) stop("Fewer than 10 eligible Reactome pathways")
  set.seed(seed); selected <- sort(sample(eligible,10))
  save_stage(list(sets=sets[selected], version=db$version), "reactome", directory)
}
if (sys.nframe()==0) reactome_stage()
