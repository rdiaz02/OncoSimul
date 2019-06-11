
## Using as model Rttf2pt1 and tth

## Copy MAGELLAN fl_* binaries to /exec/ or /exec/$R_ARCH
## We do not compile fl_draw or any of the *cgi


## binaries <- c("fl_statistics", "fl_generate", "fl_genchains")
binaries <- c("fl_statistics", "fl_generate")
oncosimul_lib <- "OncoSimulR.so"

if (WINDOWS) {
    binaries <- paste0(binaries, ".exe")
    oncosimul_lib <- "OncoSimulR.dll"
}

print(R_ARCH)
if(nzchar(R_ARCH)) {
    dest_exec <- file.path(R_PACKAGE_DIR,  paste0('exec', R_ARCH))
    dest_lib  <- file.path(R_PACKAGE_DIR,  paste0('libs', R_ARCH))
} else {
    dest_exec <- file.path(R_PACKAGE_DIR, 'exec')
    dest_lib  <- file.path(R_PACKAGE_DIR, 'libs')
}

message("Installing ", paste(binaries, collapse = " "), " to ", dest_exec)
dir.create(dest_exec, recursive = TRUE, showWarnings = FALSE)
file.copy(binaries, dest_exec, overwrite = TRUE)

message("Installing library ", oncosimul_lib, " to ", dest_lib)
dir.create(dest_lib, recursive = TRUE, showWarnings = FALSE)
file.copy(oncosimul_lib, dest_lib, overwrite = TRUE)

## Now, run file by locating it as
## system.file(package = "OncoSimulR", "exec", "fl_generate")
## or with exe under windoze


## Clean up. Otherwise, R CMD check complaints about object files
## and executables under src
try(file.remove(binaries))
try(file.remove("liblandscape.a"))
try(file.remove(list.files(file.path(getwd(), "FitnessLandscape"),
                           pattern = glob2rx("*.o"), full.names = TRUE)))
