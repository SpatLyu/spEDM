register_generic = \(name, def = NULL) {
  if (is.null(def)) {
    def = eval(bquote(function(data, ...) standardGeneric(.(name))))
  }
  methods::setGeneric(name, def)
}

for (gen in c("embedded", "fnn", "slm", "simplex", "smap",
              "multiview", "gccm", "gcmc", "sc.test")) {
  register_generic(gen)
}
