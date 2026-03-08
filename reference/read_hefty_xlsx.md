# Read time, temperature and GOF data from HeFTy output from excel file

Read time, temperature and GOF data from HeFTy output from excel file

## Usage

``` r
read_hefty_xlsx(fname)
```

## Arguments

- fname:

  path to the excel spreadsheet that contains the HeFTy outputs, i.e.
  the t-T-paths in sheet 1 and the GOF values in sheet 2

## Value

`data.frame` of the combined data

## See also

[`read_hefty()`](https://tobiste.github.io/thermoclustr/reference/read_hefty.md)

## Examples

``` r
if (FALSE) { # \dontrun{
path2myfile <- system.file("s14MM_v1.xlsx", package = "HeFTy.SmoothR")
read_hefty_xlsx(path2myfile)
} # }
```
