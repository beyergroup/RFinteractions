pdfandpng = function(expr, path, width = 5, height = 5, pointsize = 12, res = 200){
  png(paste0(path, ".png"), width = width, height = height, units = "in",
      pointsize = pointsize, res = res)
  eval(expr)
  dev.off()
  pdf(paste0(path, ".pdf"), width = width, height = height,
      pointsize = pointsize)
  eval(expr)
  dev.off()
}
