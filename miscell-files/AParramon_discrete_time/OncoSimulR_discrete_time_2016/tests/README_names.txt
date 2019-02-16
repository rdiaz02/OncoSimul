The names of the files are as self explanatory as I could think of. Some
files are names as "test.Z-somestring.R". Those are files where we set the
seed of the random number generator to a specific value.

In files with a name "test.Z*", at the very end we do a "set.seed(NULL)"
to re-initialize the seed of the random number generator. This is done so
that any files that are run after it do not have the same stream of random
numbers. Since I want to avoid re-initilizing in every file, I move those
specific test cases with fixed seeds to separate files.

Why the Z? It is true that the documentation for testthat asks that we do
not rely in a specific order for how files are tested. In this sense, any
file that uses a fixed seed cleans up after itself, so the order in which
files are tested should never mean we get a fixed stream. Nevertheless,
using the "Z" will often mean we reinitialize the seed at the very end of
all tests.
