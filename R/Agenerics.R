register_generic = \(name, def = NULL) {
  if (!methods::isGeneric(name)){
    if (is.null(def)) {
      def = eval(bquote(function(data, ...) standardGeneric(.(name))))
    }
    methods::setGeneric(name, def)
  }
}

for (gen in c("embedded", "fnn", "slm", "simplex", "smap", "ic", "multiview",
              "sc.test", "gccm", "gpc", "gcmc", "scpcm")) {
  register_generic(gen)
}
