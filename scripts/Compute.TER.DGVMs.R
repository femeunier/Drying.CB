r1 <- rast("~/Downloads/data/JULES/JULES_S2_v6.3_rh.nc")
r2 <- rast("~/Downloads/data/JULES/JULES_S2_v6.3_gpp.nc")
r3 <- rast("~/Downloads/data/JULES/JULES_S2_v6.3_npp.nc")

ter <- r1 + (r2 - r3)


writeCDF(ter, "~/Downloads/data/JULES/JULES_S2_v6.3_ter.nc", overwrite = TRUE)

#################################

r1 <- rast("~/Downloads/data/OCN/OCN_S2_ra.nc")
r2 <- rast("~/Downloads/data/OCN/OCN_S2_rh.nc")

ter <- r1 + r2

# Save result as NetCDF
writeCDF(ter, "~/Downloads/data/OCN/OCN_S2_ter.nc", overwrite = TRUE)

#################################

r1 <- rast("~/Downloads/data/ORCHIDEE/ORCHIDEE_S2Prog_1960_202412_rh.nc")
r2 <- rast("~/Downloads/data/ORCHIDEE/ORCHIDEE_S2Prog_1960_202412_ra.nc")

ter <- r1 + r2

# Save result as NetCDF
writeCDF(ter, "~/Downloads/data/ORCHIDEE/ORCHIDEE_S2Prog_1960_202412_ter.nc", overwrite = TRUE)
