.onAttach <- function(libname, pkgname) {
  packageStartupMessage("This is package ", pkgname, ".",
                        "If you are running it on an aarch64 (arm64) ",
                        "platform with a MacOS ",
                        "note that the package fails some tests ",
                        "that I have no way of debugging. ",
                        "Please read file ",
                        "README_tests_kjohnson3_aarch64-apple-darwin20.txt ",
                        "in the tests directory.")
}
