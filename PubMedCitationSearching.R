library(xml2)
library(rvest)
library(data.table)
api_key <- "e8beecd6d587f2ba1219371e358b2cf10108"
entrezUrl <- paste0("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?db=pubmed&rettype=count&api_key=", api_key)
humanGeneNames <- fread("./hsg_syns_210721.txt") #gene names file with alternate names for each gene
citationCount <- function(genes, modifiers = NULL) { #gene and other search term vectors
  citations <- matrix(0, length(genes), length(modifiers) + 1)
  rownames(citations) <- genes
  colnames(citations) <- c("Total", modifiers)
  for (i in 1:length(genes)) {
    #change given gene name to current name
    name <- humanGeneNames[match(genes[i], synonym), symbol]
    if (!is.na(name) & !(genes[i] == ""))
      searchTerm <- paste0("term=(", name, "+OR+", genes[i],")")
    else
      searchTerm <- paste0("term=(", genes[i], ")")
    #scrape PubMed for each search term
    download.file(paste0(entrezUrl, searchTerm), destfile = './search.html')
    startTime <- Sys.time()
    citations[i, 1] <- strtoi(html_text(html_nodes(read_xml('.search.html'), xpath = "//Count")))
    closeAllConnections()
    if (.1 - (Sys.time() - startTime) > 0)
      Sys.sleep(.1 - (Sys.time() - startTime))
    if (!is.null(modifiers)) {
      if (citations[i, 1] != 0) {
        for (j in 1:length(modifiers)) {
          download.file(paste0(entrezUrl, paste(searchTerm, gsub(" ", "+", modifiers[j]), sep = "+")), destfile = './search.html')
          startTime <- Sys.time()
          citations[i, j + 1] <- strtoi(html_text(html_nodes(read_xml('./search.html'), xpath = "//Count")))
          closeAllConnections()
          if (.1 - (Sys.time() - startTime) > 0)
            Sys.sleep(.1 - (Sys.time() - startTime))
        }
      } else
        citations[i, 2:(length(modifiers) + 1)] <- 0
    }
  }
  return(citations)
}
