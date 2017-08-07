rfsrc.news <- function(...) {
  newsfile <- file.path(system.file(package="_PROJECT_PACKAGE_NAME_"), "NEWS")
  file.show(newsfile)
}
