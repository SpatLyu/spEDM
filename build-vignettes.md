Due to the time-consuming computations involved in the vignettes of the *spEDM* package, 
it is necessary to pre-build the vignettes prior to package submission.

``` r
.prebuild_vignettes = \(name){
  out = paste0("vignettes/",name,".Rmd")
  inp = paste0(out,".orig")
  knitr::knit(inp,out)
}
```

-   Compile all vignettes at once

``` r
for (v in c("spEDM", "SSR", "SGC", "GCCM", "GCMC", "SCPCM")) {
  .prebuild_vignettes(v)
}
```

-   Build vignettes separately

``` r
.prebuild_vignettes("spEDM")
.prebuild_vignettes("SSR")
.prebuild_vignettes("SGC")
.prebuild_vignettes("GCCM")
.prebuild_vignettes("GCMC")
.prebuild_vignettes("SCPCM")
```
