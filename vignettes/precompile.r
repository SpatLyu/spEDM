.precompile_vignettes = \(name){
  out = paste0("vignettes/",name,".Rmd")
  inp = paste0(out,".orig")
  knitr::knit(inp,out)
}

vig = c("spEDM", "SSR", "SCT", "GCCM")

for (v in vig) {
  .precompile_vignettes(v)
}
