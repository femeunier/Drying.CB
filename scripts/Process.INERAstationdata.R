
library(readxl)

sheet_names <- excel_sheets("~/Downloads/INERA_Stations_Data.xlsx")

A <- read_excel("~/Downloads/INERA_Stations_Data.xlsx",sheet_names[2])
plot(A$Dry_bulb_temp_15h00-A$Wet_bulb_temp_15h00)

