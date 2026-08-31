rfsrc.news <- function(...) {
  newsfile <- file.path(system.file(package="randomForestSRC"), "NEWS.md")
  file.show(newsfile)
}
