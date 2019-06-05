
## Using as model Rttf2pt1 and tth

## Copy MAGELLAN fl_* binaries to /exec/ or /exec/$R_ARCH
## We do not compile fl_draw or any of the *cgi


binaries <- c("fl_statistics", "fl_generate", "fl_genchains")

if (WINDOWS) {
    binaries <- paste0(binaries, ".exe")
    execarch <- "exec"
} else {
    ## Default installation (from R extensions doc)
    print(R_ARCH)
    execarch <- if (nzchar(R_ARCH)) paste('exec', R_ARCH, sep='') else 'exec'
}

dest <- file.path(R_PACKAGE_DIR, execarch)

message("Installing ", binaries, " to ", dest)
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
file.copy(binaries, dest, overwrite = TRUE)
