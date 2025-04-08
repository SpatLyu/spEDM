.precompile_vignettes = \(name){
  out = paste0("vignettes/",name,".Rmd")
  inp = paste0(out,".orig")
  knitr::knit(inp,out)
}

.precompile_vignettes("SSR")
.precompile_vignettes("SCT")
.precompile_vignettes("GCCM")
.precompile_vignettes("GCMC")
