This is just so that I remember all the steps.

1. Install Rtools.

2. Install R. Very important that it does not install under "Program
Files" or similar, as things can breake later (e.g., using
BiocStyle). Just install under "C:\R\R-whatever", and make sure
"R-whatever" is <= 8 characters.

3. During the transition to Rtools3.3 and gcc-4.9 I've also seen issues
with compiler location. So what I do is copy the Makeconf under here to
$RHOME/etc/i386.

4. Install the required stuff in R. Remember I have issues with https (as
I am running from VirtualBox with an old Win XP). Copy the lines in
"Installs-for-OncoSimul" (including choosing the http cran mirror).


5. Of course, modify the PATH to have the new R accessible from console.

6. Copy the .tar.gz of the package to, say, C:\Documents and
settings\Ramon, and decompress.


7. Build and check as in BioC by copying the two lines of
"testing-oncosimul-windoze".
