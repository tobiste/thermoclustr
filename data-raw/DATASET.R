## code to prepare `DATASET` dataset goes here

# fname <- "inst/s14MM_v1.xlsx"

# s14MM_v1_sheet1 <- readxl::read_xlsx(fname, sheet = 1, col_names = FALSE)
# s14MM_v1_sheet2 <- readxl::read_xlsx(fname, sheet = 2, col_names = TRUE)
#
# usethis::use_data(s14MM_v1_sheet1, overwrite = TRUE)
# usethis::use_data(s14MM_v1_sheet2, overwrite = TRUE)

# tT_paths1 <- read_hefty_xlsx(fname)
# usethis::use_data(tT_paths1, overwrite = TRUE, compress = "xz")

fname1 <- 'inst/112-9-30-zr-inv.txt'
fname2 <- 'inst/112-74_50_H3_50-inv.txt'
tais1 <- read_hefty(fname1)
tais2 <- read_hefty(fname2)
saveRDS(tais1, 'inst/112-9-30-zr-inv.rds')
saveRDS(tais2, 'inst/74_50_H3_50-inv.rds')

tT_paths1 <- tais1
usethis::use_data(tT_paths1, overwrite = TRUE, compress = "xz")

tT_paths2 <- tais2
usethis::use_data(tT_paths2, overwrite = TRUE, compress = "xz")

