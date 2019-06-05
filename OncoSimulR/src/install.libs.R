
## Using as model Rttf2pt1 and tth

## Copy MAGELLAN fl_* binaries to /exec/ or /exec/$R_ARCH
## We do not compile fl_draw or any of the *cgi


binaries <- c("fl_statistics", "fl_generate") ## , "fl_genchains")

if (WINDOWS) {
    binaries <- paste0(binaries, ".exe")
    execarch <- "exec"
    libsarch <- "libs"
    oncosimul_lib <- "OncoSimulR.dll"
    ## zz: FIXME
} else {
    ## Default installation (from R extensions doc)
    print(R_ARCH)
    execarch <- if (nzchar(R_ARCH)) paste('exec', R_ARCH, sep='') else 'exec'
    libsarch <- if (nzchar(R_ARCH)) paste('libs', R_ARCH, sep='') else 'libs'
    oncosimul_lib <- "OncoSimulR.so"
}

dest <- file.path(R_PACKAGE_DIR, execarch)
dest_lib <- file.path(R_PACKAGE_DIR, libsarch)

message("Installing ", paste(binaries, collapse = " "), " to ", dest)
dir.create(dest, recursive = TRUE, showWarnings = FALSE)
file.copy(binaries, dest, overwrite = TRUE)

message("Installing library ", oncosimul_lib, " to ", dest_lib)
dir.create(dest_lib, recursive = TRUE, showWarnings = FALSE)
file.copy(oncosimul_lib, dest_lib, overwrite = TRUE)

## Now, run file by locating it as
## system.file(package = "OncoSimulR", "exec", "fl_generate")
## or with exe under windoze
