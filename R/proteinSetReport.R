

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

#' Uniprot entry source
#' 
#' Given a set of UniProt identifiers, it return the same list indicating whether the entries
#' belong to Swissprot or TrEMBL.
#' 
#' @param unimart The UniProt biomart to use for this operation
#' @param prots A vector/list of proteins to query for
#' 
#' @return A data.table containing the identifier - source database mapping (Accession,Status)
uniprotEntrySource<-function(unimart,prots) {
  return(obtainDataFromMart(mart = unimart,filters=c("accession"),values = prots, attributes = c("accession","entry_type"), key = c("accession") ))
}

#' obtain UniProt Links
#' 
#' Produces a list of uniprot identifiers for the given search field and values.
#' The search field can be ENSEMBL_Transcripts, ENSEMBL_Genes, and then the values
#' would be vectors containing either transcript identifiers or genes.
#' 
#' @param mart The biomart connection to use.
#' @param field The field to query for.
#' @param values A vector with the identifiers to query for.
#' 
#' @return A data.table containing the query to uniprot accession (and the database they belong to)
obtainUniprotLinks<-function(mart,field,values) {
  obtainDataFromMart(mart=mart,filters=c(field),
                     attributes=c(field,'uniprot_swissprot'),
                     values=values,
                     key=c(field))->uniprotSPLinks
  setnames(uniprotSPLinks,old=c('uniprot_swissprot'),new=c('uniprot'))
  uniprotSPLinks$db<-"SwissProt"
  
  uniprotTranscriptLinks<-0
  # check whether uniprotSPLinks$ensembl_transcript_id[uniprotSPLinks$uniprot_swissprot_accession==''] is not empty
  if(length(uniprotSPLinks$field[uniprotSPLinks$uniprot==''])>0) {
    obtainDataFromMart(mart=ensemblBM,filters=c(field),
                       attributes=c(field,'uniprot_sptrembl'),
                       values=unique(uniprotSPLinks$field[uniprotSPLinks$uniprot=='']),
                       key=c(field))->uniprotTREMBLLinks
    setnames(uniprotTREMBLLinks,old=c('uniprot_sptrembl'),new=c('uniprot'))
    uniprotTREMBLLinks$db<-"TrEMBL"
    rbindlist(list(uniprotSPLinks[uniprot!='',],uniprotTREMBLLinks[uniprot!='']))->uniprotTranscriptLinks
  } else {
    uniprotSPLinks[uniprot!='',]->uniprotTranscriptLinks
  }
  setkey(uniprotTranscriptLinks,uniprot) 
  return(uniprotTranscriptLinks)
}

#' chooseSilacProtein
#' 
#' Given the multiple proteins available per line on a SILAC result, choose among those with the higher amount o
#' f peptides,
#' the one that is best annotated (Swissprot vs TrEMBL)
#' 
#' @param dataTable The SILAC data set.
#' @param colForProtsAccessions The column name where the protein accessions are found. Defaults to "Protein_IDs".
#' @param colForPeptideCounts The column name where peptide counts per protein are found. Defaults to "Peptide_counts_unique".
#' @param nameForNewCol The name for the new column that will contain the chosen accession. Defaults to "Blessed_Protein_ID".
#' @param separator The token separating proteins and counts in their respective column. Defaults to ";"
#' 
#' @return Same provided dataTable with the additional column.
chooseSilacProtein<-function(dataTable,colForProtsAccessions="Protein_IDs",
                             colForPeptideCounts="Peptide_counts_unique",
                             nameForNewCol="Blessed_Protein_ID",
                             separator=";") {
  useMart(biomart = "unimart",dataset = "uniprot")->unimart
  # obtain dictionary of db source
  accessionCol<-substitute(colForProtsAccessions)
  uniprotEntrySource(unimart=unimart,prots=dataTable[,unique(unlist(strsplit(get(accessionCol),";"))),])->dict
  
  pepCol<-substitute(colForPeptideCounts)
  
  chooseProt<-function(pepCol,accessionCol,dict) {
    pepCount<-as.numeric(unlist(strsplit(pepCol,separator)))
    accession<-unlist(strsplit(accessionCol,separator))
    maxCount<-max(pepCount)
    
    selectionCount<-accession[pepCount==maxCount]
    selectionSwissProt<-selectionCount[selectionCount %in% dict[entry_type=="Swiss-Prot",accession,]]
    if(length(selectionSwissProt)>0) {
      return(selectionSwissProt[1])
    } else {
      return(selectionCount[1])
    }
  }
  dataTable[,blessedProt:=chooseProt(get(pepCol),get(accessionCol),dict),by=accessionCol]
  return(dataTable)
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
#' @param useGO Boolean, whether to use the Gene Ontology categories for the enrichment analysis.
#'        Defaults to \code{TRUE}.
#' @param usePathways Boolean, whether to use pathways controlled vocabularies: KEGG Pathways, 
#' Reactome, Reactome Interaction, Panther Pathways, BioCarta, BBID. Defaults to \code{TRUE}.
#' @param useDisease Boolean, whether to use disease related controlled vocabularies: 
#' Genetic Association DB Disease, OMIM.
#' @param useDomains Boolean, whether to use InterPro domains for the enrichment analysis. 
#' Defaults to \code{TRUE}.
#' 
#' @return the resulting enrichment analysis (functional chart) from DAVID, as a data.table. 
enrichmentAnalysisDAVID<-function(davidEmail,proteinList,listName,prot2GeneName=list(),url,
                                  useGO=TRUE,usePathways=TRUE,useDisease=TRUE,useDomains=TRUE) {
  david<-DAVIDWebService$new(email=davidEmail,url=url)
  categories<-c()
  if(useGO) {
    categories<-c("GOTERM_BP_ALL","GOTERM_MF_ALL","GOTERM_CC_ALL")
  }
  if(usePathways) {
    categories<-c(categories,"BBID","BIOCARTA","EC_NUMBER",
                  "KEGG_PATHWAY","REACTOME_PATHWAY","REACTOME_INTERACTION","PANTHER_PATHWAY")
  }
  if(useDisease) {
    categories<-c(categories,"GENETIC_ASSOCIATION_DB_DISEASE","OMIM_DISEASE")
  }
  if(useDomains) {
    categories<-c(categories,"INTERPRO")
  }
  setAnnotationCategories(david, categories)
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
  if(!is.null(foldChanges)) {
    data.frame(foldChanges=foldChanges)->foldChanges
    rownames(foldChanges)<-geneProtIdent
  }
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

#' Convert GMT Bader Lab For Piano
#'
#' Adapts the a Bader Lab (http://download.baderlab.org/EM_Genesets/) gene set file (gmt)
#' so that the gene set can be used by Piano.
#' 
#' @param pathToGMTFile The path to the .gmt file obtained from Bader lab. The file needs 
#' to be ended with .gmt extension.
#' 
#' @return the path to the converted file (same file, but .gmt extension changed to .piano.gmt)
convertBaderLabGMTForPiano<-function(pathToGMTFile) {
  readLines(pathToGMTFile)->gmtLines
  gsub("(.*?)%(.*?)%(.*?)\t(.*?)\t(.*?)","\'\\1::\\3\'\t\\2\t\\5",perl = TRUE,gmtLines)->gmtLinesMod
  gsub("\\.gmt$","\\.piano\\.gmt",pathToGMTFile)->outPath
  writeLines(gmtLinesMod,con=outPath)
  return(outPath)
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


