# Functions to run the Trans-proteomic pipeline
library(tidyverse)
library(xml2)
#library(MSGFplus) # MSGF+ currently has some bugs

print("init")

convertRaw <- function(raw.files) {
  # Convert raw files to mzML using msConvert
  
  # Template conversion command
  Convert_cmd = 'msconvert.exe MYSAMPLE.raw -o DIRNAME --filter "peakPicking true 2" --filter "msLevel 2" --zlib'
  # peakPicking true 2 centroids MS2 peaks
  # msLevel2 keeps only MS2 data (we don't need MS1 for ident or quant)
  # zlib uses compression to reduce file size
  
  # Replace the appropriate strings and run msconvert
  for (raw.file in raw.files) {
    raw.dir <- paste0(dirname(raw.file), '/')
    raw.file <- str_replace_all(raw.file, '\\\\', '/') # Fix double backslashes from choose.files()
    cmd <- str_replace(Convert_cmd, 'MYSAMPLE.raw', raw.file) %>% str_replace('DIRNAME', raw.dir)
    system(cmd)
  }
}

genCometParams <- function(project.dir, db, tol='20', termini=1, missed.cleavages=2,
                           frag.bin.tol=0.02, frag.bin.offset=0, decoy.prefix="DECOY_"){
  # Write appropriate comet search parameters
  comet.template <- "comet.params.template"
  comet.params <- paste0(project.dir,'comet.params')
  print(paste('Writing to', comet.params))
  
  cp <- readLines(comet.template)
  print(cp[1])
  cp <- gsub(pattern='\\{database\\}', replace=db, x=cp)
  cp <- gsub(pattern='\\{tolerance\\}', replace=tol, x=cp)
  cp <- gsub(pattern='\\{termini\\}', replace=termini, x=cp)
  cp <- gsub(pattern='\\{mc\\}', replace=missed.cleavages, x=cp)
  cp <- gsub(pattern='\\{frag_bin_tol\\}', replace=frag.bin.tol, x=cp)
  cp <- gsub(pattern='\\{frag_bin_offset\\}', replace=frag.bin.offset, x=cp)
  cp <- gsub(pattern='\\{decoy_prefix\\}', replace=decoy.prefix, x=cp)
  
  writeLines(cp, con=comet.params)
  print('Written')
}

run.Comet <- function(project.dir, mzML.files){
  # Run the comet search engine
  comet_cmd <- paste0("comet.exe -P", project.dir, 'comet.params ', paste(mzML.files, collapse=' '))
  system(comet_cmd)
}

# These functions "work" to run MSGF+, but there is currently a bug in the way
# MSGF+ writes out its results that breaks label-free quantification
# I need to write a subroutine that corrects the data file locations in the
# pepXML output

genMSGFParams <- function(db, tol='20 ppm', instrument=3, termini=1){
  # Write appropriate MSGF+ search parameters
  msgf.params <- msgfPar(
    database = db,
    tolerance = tol, # Parent mass tolerance
    instrument = instrument, # See instrument info below
    chargeRange = c(2,6),
    lengthRange = c(6,25),
    enzyme = 'Trypsin', # See enzyme info below
    fragmentation = 0, # 0: Auto-detect, 1: CID, 2: ETD, 3: HCD
    ntt = termini # Number of tolerable termini. 2 = fully tryptic
  )
  
  # INSTRUMENT
  # For "hybrid" spectra with high-precision MS1 and low-precision MS2, use 0.
  # For usual low-precision instruments (e.g. Thermo LTQ), use 0.
  # If MS/MS fragment ion peaks are of high-precision (e.g. tolerance = 10ppm), use 2.
  # For TOF instruments, use 2.
  # For Q-Exactive HCD spectra, use 3.
  # For other HCD spectra, use 1.
  
  # ENZYME
  # 0: unspecific cleavage, 1: Trypsin (default), 2: Chymotrypsin, 3: Lys-C, 4: Lys-N, 
  # 5: glutamyl endopeptidase (Glu-C), 6: Arg-C, 7: Asp-N, 8: alphaLP, 9: no cleavage
}

run.MSGF <- function(mzML.files, msgf.params){
  # Run the MSGF+ search engine
  # !!! Always needs a vector of mzML files, even if you only search one !!!
  # Do not import the results in to R
  res <- runMSGF(msgf.params, mzML.files, import=FALSE)
  
  # Convert mzID files to pepXML
  for (mzML.file in mzML.files) {
    bn <- tools::file_path_sans_ext(mzML.file) # Base name of the mzML file without extension
    mzid.dir <- dirname(mzML.file) # Directory for output
    mzid.file <- paste0(bn, '.mzid') # The name of the mzid file
    idconvert_cmd <- paste0('idconvert.exe ',mzid.file,' --pepXML -o ',mzid.dir)
    system(idconvert_cmd)
    unlink(mzid.file) # Delete the mzid file. We don't use it for anything else
    file.rename(paste0(bn,'.pepXML'), paste0(bn,'.pep.xml'))
  }
  
}

runTPP <- function(pepXML.files, database, combine, output.name) {
  # Run the TPP on pepXML file(s) and do label-free quant
  # Generates single protXML output
  
  # Parameterized command strings for the TPP modules
  Interact_cmd <- "InteractParser.exe MYSAMPLE.interact.pep.xml"
  Refresh_cmd <- "RefreshParser.exe MYSAMPLE.interact.pep.xml"
  Peptide_cmd <- "PeptideProphetParser.exe MYSAMPLE.interact.pep.xml ACCMASS PPM MINPROB=0 DECOYPROBS DECOY=DECOY NONPARAM EXPECTSCORE"
  iProphet_cmd <- "InterProphetParser.exe MYSAMPLE.interact.pep.xml MYSAMPLE.interact.iproph.pep.xml"
  Protein_cmd <- "ProteinProphet.exe MYSAMPLE.interact.iproph.pep.xml MYSAMPLE.interact.iproph.prot.xml NORMPROTLEN IPROPHET MINPROB=0"
  StPeter_cmd <- "StPeter.exe -d MYSAMPLE.interact.iproph.prot.xml"
  
  if(combine == TRUE){
  # Append output directory to output name
  pepxml.dir <- dirname(pepXML.files[1])
  output.name <- paste0(pepxml.dir, '/', output.name)
  
  # Run InteractParser, combining all pepXMLs
  interact_pepXML <- paste0(output.name, '.interact.pep.xml')
  interact <- paste(str_replace_all(Interact_cmd, 'MYSAMPLE', output.name), paste(pepXML.files, collapse=' '))
  system(interact)
  
  # Run RefreshParser
  refresh <- paste(str_replace(Refresh_cmd, 'MYSAMPLE', output.name), database)
  system(refresh)
  
  # Run PeptideProphet
  peptide <- str_replace(Peptide_cmd, 'MYSAMPLE', output.name)
  system(peptide)
  
  # PeptideProphet persistenly screws up the paths, which interferes
  # with StPeter. The follows is a super dirty/ugly hack to fix this.
  tx <- readLines(interact_pepXML, -1)
  tx <- gsub(pattern='C:/TPP/ShinyTPP', replace=pepxml.dir, x=tx)
  writeLines(tx, con=interact_pepXML)
  
  # Run iProphet
  ipro <- str_replace_all(iProphet_cmd, 'MYSAMPLE', output.name)
  system(ipro)
  
  # Run ProteinProphet
  protein <- str_replace_all(Protein_cmd, 'MYSAMPLE', output.name)
  system(protein)
  
  # Run StPeter
  stpeter <- str_replace(StPeter_cmd, 'MYSAMPLE', output.name)
  system(stpeter)} else {
    for(pepXML.file in pepXML.files) {
      # Get the root name of the file
      root.name <- str_replace_all(pepXML.file, '\\\\', '/') %>%
        str_replace(paste0(dirname(pepXML.file),'/'), '') %>%
        str_replace('.pep.xml', '')
      
      # Append output directory to output name
      pepxml.dir <- dirname(pepXML.file)
      output.name <- paste0(pepxml.dir, '/', root.name)
      
      # Run InteractParser, combining all pepXMLs
      interact_pepXML <- paste0(output.name, '.interact.pep.xml')
      interact <- paste(str_replace_all(Interact_cmd, 'MYSAMPLE', output.name), pepXML.file)
      system(interact)
      
      # Run RefreshParser
      refresh <- paste(str_replace(Refresh_cmd, 'MYSAMPLE', output.name), database)
      system(refresh)
      
      # Run PeptideProphet
      peptide <- str_replace(Peptide_cmd, 'MYSAMPLE', output.name)
      system(peptide)
      
      # PeptideProphet persistenly screws up the paths, which interferes
      # with StPeter. The follows is a super dirty/ugly hack to fix this.
      tx <- readLines(interact_pepXML, -1)
      tx <- gsub(pattern='C:/TPP/ShinyTPP', replace=pepxml.dir, x=tx)
      writeLines(tx, con=interact_pepXML)
      
      # Run iProphet
      ipro <- str_replace_all(iProphet_cmd, 'MYSAMPLE', output.name)
      system(ipro)
      
      # Run ProteinProphet
      protein <- str_replace_all(Protein_cmd, 'MYSAMPLE', output.name)
      system(protein)
      
      # Run StPeter
      stpeter <- str_replace(StPeter_cmd, 'MYSAMPLE', output.name)
      system(stpeter)
      
    }
  }
}

extract.results <- function(protxml.files){
  # Parse a set of protXML files to pull out relevant results
  # In "tidy" format
  out.data <- tribble(~File,~ID,~Link,~Protein.name,~Pct.coverage,~Unique.peps,~dSIN)
  print('Processing')
  
  for(protxml.file in protxml.files){
    fname <- basename(protxml.file)
    x <- read_xml(protxml.file) # Read in the results
    ns <- xml_ns_rename(xml_ns(x), d1="protXML")
    protein_groups <- xml_find_all(x, "//protXML:protein_group", ns)
    for (pg.node in protein_groups) {
      protein.node <- xml_find_first(pg.node, 'protXML:protein', ns)
      ident = xml_attr(protein.node, "protein_name")
      
      annot.node <- xml_find_first(protein.node, "protXML:annotation", ns)
      desc <- xml_attr(annot.node, 'protein_description')
      
      if(str_detect(ident, 'nxp:NX_')){
        # Need to clean up ugly nextProt description
        begin <- str_split(desc, ' \\\\Gname=')[[1]][1]
        protein.name <- str_split(begin, '\\\\Pname=')[[1]][2]
        
        # And generate nextprot URL
        url <- paste0('https://www.nextprot.org/entry/',str_replace(ident, 'nxp:', ''),'/')
        
        link_html <- paste0('<a href="',url,'">',ident,'</a>')
        
      } else {
        protein.name <- desc
        link_html <- ''
      }
      coverage = as.numeric(xml_attr(protein.node, "percent_coverage"))
      distinct.peps = as.numeric(xml_attr(protein.node, "total_number_distinct_peptides"))
      stpeter.node <- xml_find_first(protein.node, "protXML:analysis_result/protXML:StPeterQuant", ns)
      dSIn <- as.numeric(xml_attr(stpeter.node, "dSIn"))
      out.row <- tribble(~File, ~ID,~Link,~Protein.name,~Pct.coverage,~Unique.peps,~dSIN,
                         fname, ident, link_html, protein.name, coverage, distinct.peps, dSIn)
      out.data <- bind_rows(out.data, out.row)
    }
  }
  out.data %>% filter(!is.na(dSIN)) %>% arrange(-dSIN) -> out.data
  outfile <- paste0(dirname(protxml.files[1]),'/result_table.csv')
  write_csv(out.data, path = outfile)
  print(paste('Results written to', outfile))
  return(out.data)
}

## Example command-line workflow
# 1. Get all raw files in a directory and convert them
#   Set this to the directory where your raw files are:
#> project.dir <- "E:/2016-08_StPeter_Validation/"

#   Don't edit these lines:
#> raw.files <- list.files(project.dir, pattern = 'raw', full.names=T)
#> convertRaw(raw.files)

# 2. Generate comet parameters
#   Set this to your desired search database
#> database <- "E:/dbase/nextprot_all_DECOY.fasta"

#> genCometParams(project.dir, database, frag.bin.tol=1.0005, frag.bin.offset=0.4)
# Note: if you are using high-res MS2 (e.g. HCD in QEHF),
# remove the frag.bin.tol and frag.bin.offset arguments,
# e.g. genCometParams(project.dir, database)

# 3. Search all mzML files in the project directory with Comet
#   You should not need to edit these lines
#> mzML.files <- list.files(project.dir, pattern = 'mzML', full.names=T)
#> run.Comet(project.dir, mzML.files)

# 4a. Run the TPP, COMBINING ALL SEARCH RESULTS INTO ONE OUTPUT
# This is useful for, e.g., gel band data where you need to combine many files
#> pepxml.files <- list.files(project.dir, pattern = 'pep.xml', full.names=T)
#> runTPP(pepxml.files, database, 'combined')
# Note: 'combined' can be changed to whatever name you like,
# just make sure it's in quotes. The output will be
# 'combined.interact.iproph.prot.xml'

# 4b. Run the TPP but process each search result individually
#> pepxml.files <- list.files(project.dir, pattern = 'pep.xml', full.names=T)
#> for(pepxml.file in pepxml.files){
#>  bn <- tools::file_path_sans_ext(mzML.file)
#>  runTPP(c(pepxml.file), database, bn)
#>}