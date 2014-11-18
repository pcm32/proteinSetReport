

#' Obtain data from mart
#' 
#' \code{obtainDataFromMart} encapsulates a biomart call to produce a data.table object with the result.
#' 
#' @param mart the biomart object previously initialized.
#' @param filters the set of filters to be used on the biomart call.
#' @param attributes the set of attributes to be used on the biomart call.
#' @param values the query values to be used, in equivalent dimension as the filters.
#' @param key the column to be used as key in the data.table, this should be part of the attributes
#' 
#' @return The biomart result as a data.table indexed by key.
obtainDataFromMart<-function(mart,filters,attributes,values,key) {
  return(data.table(getBM(attributes=attributes,filters=filters,values=values,mart=mart,uniqueRows=T),key=key))
}

#' Enrichment analysis DAVID
#' 
#' \code{enrichmentAnalysisDAVID} uses the \code{RDAVIDWebService} package to execute an
#' enrichment analysis over http with the DAVID server. 
#' 
#' @param davidEmail The email registered at DAVID to be able to use their web service.
#' @param proteinList The list of proteins (UniProt Accessions) that should be analyzed.
#' @param listName A label for the list.
#' @param prot2GeneName An optional dictionary (in a data.table) with prot and GeneName columns,
#'        used to translate the uniprot accessions to Gene names.
#' @return the resulting enrichment analysis (functional chart) from DAVID, as a data.table.
enrichmentAnalysisDAVID<-function(davidEmail,proteinList,listName,prot2GeneName=list()) {
  david<-DAVIDWebService$new(email=davidEmail)
  setAnnotationCategories(david, c("GOTERM_BP_ALL","GOTERM_MF_ALL","GOTERM_CC_ALL",
                                   "BBID","BIOCARTA","EC_NUMBER",
                                   "KEGG_PATHWAY","REACTOME_PATHWAY","REACTOME_INTERACTION","PANTHER_PATHWAY",
                                   "GENETIC_ASSOCIATION_DB_DISEASE","OMIM_DISEASE","INTERPRO"
  ))
  result<-addList(david,inputIds=proteinList,listName=listName,idType="UNIPROT_ACCESSION")
  getFunctionalAnnotationChart(david)->f
  dataInF<-matrix(nrow = nrow(f), ncol = ncol(f))
  for(i in 1:length(f@names)) {
    dataInF[,i]<-unlist(f@.Data[i])
  }
  names(dataInF)<-f@names
  data.table( dataInF )->dataInF.dt
  categoryNames<-categories(f)  
  setnames(dataInF.dt,f@names)
  # We add the categories names instead of the per session category indexes.
  dataInF.dt[,Category:=categoryNames[as.numeric(as.character(Category))],][,PValue:=as.numeric(as.character(PValue)),]
  dataInF.dt[,Term:=as.character(Term),][,Genes:=as.character(Genes),][,FDR:=as.numeric(as.character(FDR)),]
  
  if(length(prot2GeneName)>0) {
    # add GeneNames if the dictionary is provided.
    dataInF.dt[,GeneNames:=paste( prot2GeneName[prot %in% unlist(strsplit(as.character(Genes),split = ", ")),]$GeneName ,collapse=", "),by=Term]
  }
  return(dataInF.dt)
}

#' Get ENTREZ GI From ENSEMBL IDs
#' 
#' This function uses a biomart connection to ENSEMBL biomart to fetch ENSEMBL gene ids to
#' ENTREZ GI.
#' 
#' @param mart The biomart (ENSEMBL) instance to use
#' @param ensemblIDs A vector of ENSEMBL IDs to get ENTREZ GIs
#' 
#' @return A data.table view of the biomart result, with \code{ensembl_gene_id} as key and \code{entrezgene}.
getENTREZGIFromENSEMBLIDs<-function(mart,ensemblIDs) {
  field<-'ensembl_gene_id'
  obtainDataFromMart(mart=mart,filters=c(field),
                     attributes=c(field,'entrezgene'),
                     values=ensemblIDs,
                     key=c(field))
}

#' Get ENTREZ GI From ENSEMBL IDs
#' 
#' This function uses a biomart connection to UniProt to fetch proteins to ensembl gene ids, and a
#' ENSEMBL biomart to fetch ENSEMBL gene ids to ENTREZ GI.
#' 
#' @param martUniprot The uniprot biomart instance to use
#' @param martEnsembl The ensembl biomart instance to use
#' @param ensemblIDs A vector of UniProt IDs to get ENTREZ GIs
#' 
#' @return A data.table view of the biomart result, with \code{accession} and \code{entrezgene}.
getENTREZGIFromUNIPROTIDs<-function(martUniprot,martEnsembl,uniprotIDs) {
  field<-'accession'
  uniprot2Ensembl<-obtainDataFromMart(mart=martUniprot,filters=c(field),
                      values=uniprotIDs,attributes=c(field,'ensembl_id'),key=c('ensembl_id'))
  
  ens2gi<-getENTREZGIFromENSEMBLIDs(mart = martEnsembl, ensemblIDs = unique(uniprot2Ensembl$ensembl_id))
  
  ens2gi[uniprot2Ensembl]
}

#' Enrichment Analysis REACTOME
#' 
#' Runs enrichment analysis using the ReactomePA.
#' 
#' @param giList A vector of entrez gene identifiers to run the analysis for
#' @param gi2GeneName An optional dictionary (as a data.table) with columns gi \code{gi} to \code{GeneName}, 
#' which maps entrez gene ids to gene names.
#' 
#' @return The Reactome enrichment analysis for the vector of gis given, as a data.table
enrichmentAnalysisREACTOME<-function(giList,gi2GeneName=list()) {
  enrichPathway(gene = giList, organism = "human", pAdjustMethod = "fdr", readable = F, pvalueCutoff = 0.9, qvalueCutoff = 0.9)->rpa
  data.table(rpa@result)->res
  if(length(gi2GeneName)>0) {
    res[,GeneNames:=paste( gi2GeneName[gi %in% unlist(strsplit(as.character(geneID),split = "/")),]$GeneName, collapse=", "),by=ID]
  }
}

#' Parse DAVID KEGG Pathway Term
#' 
#' Parses a "term" from the result of a DAVID enrichment analysis, as produced by function \code{enrichmentAnalysisDAVID},
#' to extract KEGG Pathways enriched to be feed to the pathview function to colour KEGG Pathway diagrams.
#' 
#' @param term An individual term from the DAVID enrichment analysis. This is not a vectorized function.
#' 
#' @return A named list with species, pathID and pathName. These are the required fields for the pathview function, 
#' besides the list of proteins to colour.
parseDavidKEGGPathwayTerm<-function(term) {
  # parses this kind of input hsa03010:Ribosome
  str_match_all(string = term, pattern = "(\\D+)(\\d+):(.*)")->parsed
  if(!is.na(parsed) && length(parsed) == 1 && length(parsed[[1]]==4)) {
    return(list(species=parsed[[1]][1,2],pathID=parsed[[1]][1,3],pathName=parsed[[1]][1,4]))
  }
  return(list(species=NA,pathID=NA,pathName=NA))
}

#' Make link for DB
#' 
#' Lower level function of the link making series. 
#' 
#' @param entry Identifier of the database where we are linking to
#' @param urlPrefix The prefix to be used before the entry for constructing the URL.
#' @param urlPostfix The postfix to be used after the entry for constructing the URL. Defaults to an empty string.
#' @param dict An optional data.frame which has two columns \code{entry} and \code{name}, and acts as a dictionary to
#' translate the identifier in entry to something else (such as a Gene name) to be used as label of the link.
#' @param type Either "html" or "markdown", defining the type of link generated.
#' 
#' @return the string with the link.
#' 
makeLinkForDB<-function(entry,urlPrefix,urlPostfix="",dict,type) {
  if(is.null(entry) || is.na(entry) || str_length(entry)==0 ) {
    return("")
  }
  label<-entry
  if(is.data.frame(dict)) {
    if(length(entry)>1) {
      label<-dict$name[dict$entry %in% entry]
    } else {
      label<-dict$name[dict$entry==entry] 
    }
  }
  if(type=="html") {    
    return(paste("<a href=\"",urlPrefix,unlist(lapply(entry,FUN = URLencode, reserved=T)),urlPostfix,"\">",label,"</a>",sep=""))
  } else if(type=="markdown") {
    return(paste("[",label,"](",urlPrefix,unlist(lapply(entry,FUN = URLencode, reserved=T)),")",sep=""))
  }
}

#' Make UniProt link
#' 
#' Produces a link to UniProt for the specified entry.
#' 
#' @param uniprotEntry The accession of the UniProt entry.
#' @param dict The dictionary (data frame with entry and name columns) to 
#' translate the entry into something different for the URL label. Defaults to \code{NA}.
#' @param type Either "html" or "markdown", defining the type of link generated. Defaults to "html".
#' 
#' @return The link for the UniProt entry.
#' 
#' @seealso \code{\link{makeLinkForDB}}
#' 
makeUniprotLink<-function(uniprotEntry,dict=NA,type="html") {
  prefix<-"http://www.uniprot.org/uniprot/"
  makeLinkForDB(entry=uniprotEntry,urlPrefix=prefix,dict=dict,type=type)
}

#' Make REACTOME link
#' 
#' Produces a link to Reactome for the specified entry.
#' 
#' @param reactomeEntry The accession of the Reactome pathway entry.
#' @param dict The dictionary (data frame with entry and name columns) to 
#' translate the entry into something different for the URL label. Defaults to \code{NA}.
#' @param type Either "html" or "markdown", defining the type of link generated. Defaults to "html".
#' 
#' @return The link for the Reactome Pathway entry.
#' 
#' @seealso \code{\link{makeLinkForDB}}
#' 
makeReactomeLink<-function(reactomeEntry,dict=NA,type="html") {
  prefix<-"http://www.reactome.org/cgi-bin/link?SOURCE=Reactome&ID="
  makeLinkForDB(entry=reactomeEntry,urlPrefix=prefix,dict=dict,type=type)
}

#' Make InterPro protein architecture (UniProt) link
#' 
#' Produces a link to InterPro for a UniProt specified entry.
#' 
#' @param uniprotEntry The accession of the UniProt entry.
#' @param type Either "html" or "markdown", defining the type of link generated. Defaults to "html".
#' 
#' @return The link to the InterPro page for the UniProt entry.
#' 
#' @seealso \code{\link{makeLinkForDB}}
#'
makeInterProtDomainArchLink<-function(uniprotEntry,type="html") {
  prefix<-"http://www.ebi.ac.uk/interpro/protein/"
  makeLinkForDB(entry=uniprotEntry,urlPrefix = prefix,type=type,dict=NA )
}

#' Make InterPro link
#' 
#' Produces a link to InterPro for the specified entry.
#' 
#' @param interproEntry The accession of the UniProt entry.
#' @param dict The dictionary (data frame with entry and name columns) to 
#' translate the entry into something different for the URL label. Defaults to \code{NA}.
#' @param type Either "html" or "markdown", defining the type of link generated. Defaults to "html".
#' 
#' @return The link for the InterPro domain entry.
#' 
#' @seealso \code{\link{makeLinkForDB}}
#' 
makeInterproLink<-function(interproEntry,dict=NA,type="html") {
  prefix<-"http://www.ebi.ac.uk/interpro/entry/"
  makeLinkForDB(entry=interproEntry,urlPrefix = prefix,type=type,dict=dict )
}

#' Make Molecular Interaction Ontology Link
#' 
#' \code{makeMolIntOntologyLink} Produces a link to the Ontology Look-up Service (OLS) from the EBI,
#' for the specified entry of the MI Ontology.
#' 
#' This function should probably be improved to receive a dictionary instead of the description in the
#' input.
#' 
#' @param line A line containing the accession of the MI Ontology followed by the description.
#' @param dict The dictionary (data frame with entry and name columns) to 
#' translate the entry into something different for the URL label. Defaults to \code{NA}.
#' @param type Either "html" or "markdown", defining the type of link generated. Defaults to "html".
#' 
#' @return The link for the MI.
#' 
#' @seealso \code{\link{makeLinkForDB}}
#' 
makeMolIntOntologyLink<-function(line,type="html") {
  # line is like "MI:0064"(interologs mapping)
  gsub("\"","",line)->line
  str_match_all(string = line,pattern = "(MI:\\d{4})\\((.*)\\)")->parsed
  
  checkParse<-function(parsedPart,line="none") {
    if(length(parsedPart)==3) {
      miID<-parsedPart[1,2]
      desc<-parsedPart[1,3]
      data.frame(name = c(desc), entry = c(miID))->dict
      prefix<-"http://www.ebi.ac.uk/ontology-lookup/browse.do?ontName=MI&termId="
      makeLinkForDB(entry = miID, urlPrefix = prefix, dict = dict, type = type)
    } else {
      line
    }
  }
  
  lapply(parsed,FUN = checkParse, line = line )
}


