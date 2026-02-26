## code to prepare `DATASET` dataset goes here

# fname <- "inst/s14MM_v1.xlsx"

# s14MM_v1_sheet1 <- readxl::read_xlsx(fname, sheet = 1, col_names = FALSE)
# s14MM_v1_sheet2 <- readxl::read_xlsx(fname, sheet = 2, col_names = TRUE)
#
# usethis::use_data(s14MM_v1_sheet1, overwrite = TRUE)
# usethis::use_data(s14MM_v1_sheet2, overwrite = TRUE)

# tT_paths1 <- read_hefty_xlsx(fname)
# usethis::use_data(tT_paths1, overwrite = TRUE, compress = "xz")

fname1 <- "inst/112-73_30_H1_50-inv.txt"
# fname2 <- "inst/112-72_30_H1_50-inv.txt"
tais1 <- read_hefty(fname1)
# tais2 <- read_hefty(fname2)
saveRDS(tais1, "inst/112-73.rds")
# saveRDS(tais2, "inst/112-72_30_H1_50-inv.rds")

all(validUTF8(tT_paths1$summary$grain))
# all(validUTF8(tT_paths2$summary$grain))
tT_paths1$summary$grain <- iconv(tT_paths1$summary$grain, from = "", to = "ASCII//TRANSLIT", sub = "")
# tT_paths2$summary$grain <- iconv(tT_paths2$summary$grain, from = "", to = "ASCII//TRANSLIT", sub = "")
all(validUTF8(tT_paths1$summary$grain))
# all(validUTF8(tT_paths2$summary$grain))

tT_paths <- tais1 |> crop_paths(time = c(0, 400), temperature = c(0, 250))
usethis::use_data(tT_paths, overwrite = TRUE)

# tT_paths2 <- tais2 |> crop_paths(time = c(0, 400), temperature = c(0, 250))
# usethis::use_data(tT_paths2, overwrite = TRUE)
