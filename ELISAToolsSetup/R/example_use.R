
# Import library ----------------------------------------------------------
library(ELISAToolsSetup)



# Clear Environment -------------------------------------------------------
rm(list = ls())





# Pipeline Wrapper -------------------------------------------------------






# Test -----------------------------------------------------------------
add.plate.replicate.IDs(plate = plate)
add.plate.replicate.IDs(plate = plate.april)

mfi %>% create.id_replicate(id.type = "sample")
mfi.april %>% create.id_replicate(id.type = "specimen")



# Extract plates MFI only -----------------------------------------------------
#Craig plates
temp <- mfi %>% extract.plate() %>%
  create.id_replicate(id.type = "sample") %>%
  loop.mfi.sub(plate = add.plate.replicate.IDs(plate = plate))

#April plates
temp.april <- mfi.april %>% filter(plate == 1) %>%
  create.id_replicate(id.type = "specimen") %>%
  loop.mfi.sub(plate = add.plate.replicate.IDs(plate = plate.april))


# Extract OD formatted plates ---------------------------------------------
#Craig full pipeline
temp.2 <- mfi %>% extract.plate() %>%
  create.id_replicate(id.type = "sample") %>%
  loop.mfi.sub(plate = add.plate.replicate.IDs(plate = plate)) %>%
  format.analyte.plates(plate.batch = 3)

#April full pipeline
temp.april.2 <- mfi.april %>% filter(plate == 1) %>%
  create.id_replicate(id.type = "specimen") %>%
  loop.mfi.sub(plate = add.plate.replicate.IDs(plate = plate.april)) %>%
  format.analyte.plates(plate.batch = 1)


# Export ------------------------------------------------------------------
#Export Craig
OD.export(temp.2,
          export.dir = "/Users/sar/Dropbox/Wellness analysis 2022/ODformat_plates/craig_plate3",
          plate.batch = 3)

#Export April
OD.export(temp.april.2,
          export.dir = "/Users/sar/Dropbox/Wellness analysis 2022/ODformat_plates/april_plate1_th",
          plate.batch = 1)
