# Load libraries ----------------------------------------------------------
library(tidyverse)


# Pipeline Sub-Functions -------------------------------------------------
# === Wrangle Plate Data ===
add.replicate.id <- function(x, #value
                             y #value in previous column (-1)
){
  if(x == y){
    new.replicate.number <- (str_split(x, "_")[[1]][2] %>% as.numeric) + 1
    out <- gsub(pattern = "_([0-9])", replacement = paste0("_", new.replicate.number), x = x)
  } else
  {
    out <- x
  }

  return(out)
}

add.1 <- function(x){
  out <- paste0(x, "_1")
  return(out)
}

remove.whitespace <- function(x){
  out <- gsub(pattern = " ", replacement = "", x = x)
  return(out)
}

loop.replicates <- function(plate.df,
                            orientation #one of h or v; declare whether replicates have h = horizontal, v = vertical orientation
){
  if(!(orientation %in% c("h", "v"))){
    stop("Invalid orientation input")
  }

  if(orientation == "h"){
    out <- data.frame(V1 = plate.df[,1])
    for(i in 2:ncol(plate.df)){
      out[,i] <- mapply(add.replicate.id, x = plate.df[,i],y = plate.df[,i - 1])
    }
  }



  if(orientation == "v"){
    stop("Vertical orientation of replicates is not supported in the current version.")
  }


  return(out)
}

extract.plate <- function(mfi, well.name = well){
  out <- mfi %>% mutate(plate = str_extract({{well.name}},"P\\d+")) %>% #extract P and any number of following digits
    mutate(plate = gsub(pattern = "P", replacement = "", x = plate))
  return(out)
}

force.bind <- function(df1, df2) {
  colnames(df2) = colnames(df1)
  bind_rows(df1, df2)
}

data.frame.OD.rows <- function(row){
  out <- row %>% data.frame %>% t() %>% data.frame
  return(out)
}

# === Substituting specimen/sample ID with MFI ===
extract.specimen <- function(mfi, well.name = well){
  out <- mfi %>% mutate(specimen = sapply(str_split({{well.name}}, "_"), "[[", 3) %>%
                          gsub(pattern = "-([0-9])", replacement = "")
  )
  return(out)
}


# Pipeline Functions ------------------------------------------------------
#1. Add replicate IDs to plate
add.plate.replicate.IDs <- function(plate){
  out <- plate %>%
    mutate(across(everything(), remove.whitespace)) %>% #Remove whitespace from sample names
    mutate(across(everything(), add.1)) %>% #prepare; add 1 to all samples to initialize replicate ID loop
    loop.replicates(orientation = "h") #Add replicate IDs to plate according to replicate orientation
  return(out)
}

#2. Create id_replicate column
create.id_replicate <- function(mfi,
                                id.type, #one of: "specimen" or "sample"
                                well.name = well #name of column containing .fcs file name with specimen ID embedded (i.e. Specimen_001_C0-1_P3_015.fcs). Only required if id.type == "specimen"
){
  if(!(id.type %in% c("specimen", "sample"))){
    stop("Invalid id.type provided. Must be one of: 'specimen' or 'sample'")
  }

  if(id.type == "specimen"){
    out <- mfi %>% extract.specimen() %>% mutate(id_replicate = paste(specimen, replicate, sep = "_"))
  }

  if(id.type == "sample"){
    out <- mfi %>% mutate(id_replicate = paste(sample, replicate, sep = "_"))
  }


  return(out)
}


#3. Susbstitute specimen/sample ID with MFI on plate map
loop.mfi.sub <- function(mfi, #with id_replicate column (run #1: create.id_replicate)
                         plate, #with wells containing IDs and replicate (run #2: add.plate.replicate.IDs)
                         metadata = c("experiment", "well", "sample", "dilution", "replicate", "sample_type", "specimen", "plate") #all columns except analytes
)
{
  # === identify analyte names ===
  analytes <- mfi %>% select(-c(any_of(metadata), id_replicate)) %>% colnames

  # === subsitute speciment/sample ID with MFI
  #Initialize loop variables
  out <- list()
  well_row <- c("A", "B", "C", "D", "E", "F", "G", "H")

  #Loop
  for(i in 1:length(analytes)){
    analyte <- analytes[i]
    print(paste("Creating plate for:", analyte))

    analyte.plate.df <- data.frame(rows = rownames(plate))
    mfi.analyte <- mfi %>% select(id_replicate, analyte)
    for(j in 1:ncol(plate)){
      plate.column <- data.frame(id_replicate = plate[,j]) %>% left_join(mfi.analyte, by = "id_replicate")
      analyte.plate.df[,j+1] <- plate.column[,2]
    }

    #clean-up plate df
    analyte.plate.df <- analyte.plate.df %>% select(-rows)
    colnames(analyte.plate.df) <- colnames(plate)

    #re-format plate in long-form  (col 1 = well, col 2 = mfi)
    analyte.plate.longer <- analyte.plate.df %>% cbind(well_row) %>%
      pivot_longer(cols = colnames(analyte.plate.df), names_to = "well_col", values_to = "mfi") %>%
      mutate(well_col = gsub("V", "", well_col)) %>% mutate(well = paste0(well_row, well_col)) %>%
      select(well, mfi)

    #save long-form analyte plate to list
    out[[i]] <- analyte.plate.df
  }
  names(out) <- analytes



  return(out)
}


#4. Format as OD files (version 1)
format.analyte.plates <- function(analytes.list,
                                  plate.batch #plate number (aka batch number)
)
{
  # === constant rows between plates ===
  row.1 <- c("##BLOCKS= 1", rep("", 13)) %>% data.frame.OD.rows
  row.3 <- c("", "Temperature(iC", 1:12) %>% data.frame.OD.rows
  row.12 <- c("~End", rep("", 13)) %>% data.frame.OD.rows


  # === Add rows to each plate ===
  #initialize loop
  out <- list()
  for(i in 1:length(analytes.list)){
    analyte <- names(analytes.list[i])
    print(paste("OD formatting", analyte))

    analyte.df <- analytes.list[i] %>% data.frame

    #rows
    row.2 <- c("Plate:", paste(plate.batch, analyte), "PlateFormat", "Endpoint", "Absorbance", "Raw", FALSE, 1, rep("", 6)) %>% data.frame.OD.rows
    row.4 <- cbind("", "24", analyte.df[1,]) %>% data.frame
    row.5_11 <- cbind(rep("", nrow(analyte.df[-1,])),
                      rep("", nrow(analyte.df[-1,])),
                      analyte.df[-1,]) %>% data.frame

    #Loop rbind each row
    row_names <- c("row.1", "row.2", "row.3", "row.4", "row.5_11", "row.12", "row.13", "row.14")
    #initialize loop
    analyte.df.out <-row.1
    for(k in 2:length(row_names)){
      analyte.df.out <- force.bind(analyte.df.out,
                                   get(row_names[k]) %>% mutate_all(as.character))
    }
    rownames(analyte.df.out) <- NULL
    out[[i]] <- analyte.df.out

  }





  names(out) <- names(analytes.list)
  return(out)
}



#5. Export files
OD.export <- function(ODplates, #output from format.analyte.plates
                      export.dir, #filepath
                      plate.batch)
                      {
  for(i in 1:length(ODplates)){

    #prepare
    export.df <- ODplates[i]
    export.filename <- paste0(names(ODplates[i]), "_plate", plate.batch, ".txt")
    export.full.path <- paste(export.dir, export.filename, sep = "/")

    #export
    write.table(x = export.df, file = export.full.path,
                sep = "\t",
                quote = FALSE, col.names = FALSE, row.names = FALSE)
      }
}
