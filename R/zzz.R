.onAttach <- function(libname, pkgname) {
  rfsrc.version <- read.dcf(file=system.file("DESCRIPTION", package=pkgname), 
                            fields="Version")
  packageStartupMessage(paste("\n",
                              pkgname,
                              rfsrc.version,
                              "\n",
                              "\n",
                              "Type rfsrc.news() to see new features, changes, and bug fixes.",
                              "\n",
                              "\n"))
}
