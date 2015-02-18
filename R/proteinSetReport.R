

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
  
  ens2gi[uniprot2Ensembl,allow.cartesian=TRUE]
}

#' Get UniProt Protein Names
#' 
#' Retrieves the UniProt Gene name for the provided ENSEMBL identifiers. The ENSEMBL
#' Biomart data set is choosen through the species parameter. Multiple names are collapsed
#' into a comma separated list, so that ENSEMBL gene ids are unique.
#' 
#' @param ensemblIDs The identifiers to search names for. These should belong to the
#' data set specifief by \code{species}
#' @param species The name of the ENSEMBL Biomart data set to use
#' 
#' @return A data.table with rows \code{ensembl_gene_id} and \code{uniprot_genename}. 
getUniprotProteinNames<-function(ensemblIDs,species="mmusculus_gene_ensembl") {
  bm <- useMart("ensembl")
  bm <- useDataset(species, mart=bm)
  # Get ensembl gene ids and GO terms
  eg2Name <- data.table(getBM(mart=bm, attributes=c('ensembl_gene_id','uniprot_genename'), filters = c('ensembl_gene_id'), values = ensemblIDs))
  eg2Name[,uniprot_genename:=gsub(";","",uniprot_genename),]
  eg2Name[,list(name=paste(unique(uniprot_genename),collapse=", ")),by=ensembl_gene_id]->CvsT_names
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
  # deals with certain descriptions carrying non UTF-8 characters.
  res[,Description:=iconv(Description, "UTF-8", "UTF-8",sub=' '),by=ID]
  if(length(gi2GeneName)>0) {
    res[,GeneNames:=paste( gi2GeneName[gi %in% unlist(strsplit(as.character(geneID),split = "/")),]$GeneName, collapse=", "),by=ID]
  }
  return(res)
}

#' Enrichment Analysis REACTOME WS
#' 
#' Runs enrichment analysis using the REACTOME REST web service.
#' 
#' @param proteins A vector of UniProt identifiers to run the analysis for
#' 
#' @return The Reactome enrichment analysis for the vector of proteins given, as a data.table
enrichmentAnalysisREACTOMEWS<-function(proteins) {
  REACTOMEEnrichmentAnalysis(identifiers = proteins)->res
  data.table(res,key='Pathway.identifier')->res
  setnames(res,old=c('Pathway.identifier','Pathway.name','Entities.pValue','Entities.FDR','Submitted.entities.found'),
           new=c('ID','Description','pvalue','qvalue','geneID'))
  res[,ID:=as.character(ID),]
  return(res)
} 

#' Choose Enrichment Rows
#' 
#' Processes an enrichment analysis data.table (which contains a pvalue and qvalue rows)
#' to produce a selection of the rows, considering that a very high number of rows would
#' be problematic for most people in an HTML report.
#' 
#' @param dataTable A data.table object holding results from enrichment, such as from
#' \code{enrichmentAnalysisREACTOMEWS} or \code{enrichmentAnalysisDAVID}.
#' @param maxDefault The maximum number of rows to use if there are many significant results.
#' 
#' @return the shortened version of the provided dataTable.
chooseEnrichmentRows<-function(dataTable,maxDefault=30) {
  if(nrow(dataTable[pvalue<0.05,])>maxDefault) {
    if(nrow(dataTable[qvalue<0.05,])>maxDefault) {
      dataTable[1:maxDefault,]->dataTable
    } else {
      dataTable[qvalue<0.05,]->dataTable
    }
  } else {
    dataTable[1:min(maxDefault,nrow(dataTable)),]->dataTable
  }
  return(dataTable)
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

#' Run piano GSEA analysis
#' 
#' Runs the piano gene set enrichment analysis for a set of proteins and accompanying p-values.
#' 
#' @param geneProtIdent A list of proteins or gene identifiers, which matches the given gene set collection.
#' @param pvalues A list of p-values, of the same lenght as the protein list.
#' @param foldChanges 
#' @param pianoGSC A piano formatted gene set collection
#' 
#' 
#' @return A piano GSEA result object.
runPianoGSEAnalysis<-function(geneProtIdent,pvalues,foldChanges=NULL,pianoGSC,minGSSize=5,maxGSSIze=300) {
  data.frame(pvalues=pvalues)->pvaluesDF
  rownames(pvaluesDF)<-geneProtIdent
  runGSA(geneLevelStats = pvaluesDF,directions = foldChanges, gsc=pianoGSC,gsSizeLim=c(minGSSize,maxGSSIze))
}

#' Merge Genes to GSA Res
#' 
#' Piano GSA result summary table doesn't include the detail of which genes in the query are
#' part of each enriched category/group. This method joins the summary to the genes that belong
#' to each category. It can cut the list if a cut-off for the adjusted p-value is given.
#' 
#' @param gsaRes A piano GSA result object.
#' @param padjCutoff A cut-off for rows to show in the final table; defaults to \code{0.5}.
#' @param species The ensembl biomart data set name for retrieving names for the genes. This will only work
#' if the GSA piano method was executed with a data set where the identifier for genes are ENSEMBL IDs.
mergeGenesToGSARes<-function(gsaRes,padjCutoff=0.5,species=NA) {
  data.table(GSAsummaryTable(gsaRes))->tb.dt
  data.table(category=names(gsaRes$gsc),
             geneIdents=sapply(gsaRes$gsc,function(x) paste(x,collapse=", ")))->category2genes
  setkey(tb.dt,Name)
  setkey(category2genes,category)
  tb.dt[category2genes]->mergedWithGenes
  padjCutoff<-0.5
  
  selection<-rep(FALSE,times = nrow(mergedWithGenes))
  
  cols<-c('p adj (dist.dir.up)', 'p adj (dist.dir.dn)', 'p adj (non-dir.)', 'p adj (mix.dir.up)', 'p adj (mix.dir.dn)')
  
  for(col in cols) {
    if(col %in% colnames(mergedWithGenes)) {
      selection<-(selection | mergedWithGenes[[col]] <= padjCutoff)
    }
  }

  
  if(!is.na(species)) {
    unique(unlist(strsplit(mergedWithGenes$geneIdents[selection],split = ", ")))->uniqueENSEMBLIDs
    getUniprotProteinNames(uniqueENSEMBLIDs,species)->ensemblIds2GeneNames
  
    mergedWithGenes[,GeneNames:=paste(
      ensemblIds2GeneNames[ensembl_gene_id %in% unlist(strsplit(geneIdents,split = ", "))]$name,
      collapse = ", "),by=Name]
  }
  mergedWithGenes[order(as.numeric(mergedWithGenes$'p adj (non-dir.)')),]
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
      label<-dict$name[match(entry,dict$entry)]
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


