``` r
.precompile_vignettes = \(name){
  out = paste0("vignettes/",name,".Rmd")
  inp = paste0(out,".orig")
  knitr::knit(inp,out)
}
```

-   Compile all vignettes at once

``` r
for (v in c("spEDM", "SSR", "SCT", "GCCM", "GCMC", "SCPCM")) {
  .precompile_vignettes(v)
}
```

-   Build vignettes separately

``` r
.precompile_vignettes("spEDM")
.precompile_vignettes("SSR")
.precompile_vignettes("SCT")
.precompile_vignettes("GCCM")
.precompile_vignettes("GCMC")
.precompile_vignettes("SCPCM")
```
