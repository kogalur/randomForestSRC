rfsrc.news <- function(...) {
  newsfile <- file.path(system.file(package="randomForestSRC"), "NEWS")
  file.show(newsfile)
}
