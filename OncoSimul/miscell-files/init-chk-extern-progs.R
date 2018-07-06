## Testing CBN programs available
.._OncoSimul_test.ctcbn <- system("ct-cbn -h", ignore.stdout = TRUE)
.._OncoSimul_test.hcbn <- system("h-cbn -h", ignore.stdout = TRUE)
if(.._OncoSimul_test.ctcbn || .._OncoSimul_test.hcbn) {
    warning(paste(
        "\n\n",
        "\n******************************************************",
        "\n******************************************************\n",
        "          OncoSimulR installation warning:\n",
        "\n",
        "The external programs h-cbn and ct-cbn were not found.",
        "You will not be able to use CBN.",
        "You can download them from http://www.bsse.ethz.ch/cbg/software/ct-cbn",
        "You might want to use the makefile Makefile-cbn-modified-RDU under",
        "the miscell-code directory of this package.\n",
        "\n******************************************************",
        "\n******************************************************",
        "\n\n"
        )
            )
}
rm(.._OncoSimul_test.ctcbn)
rm(.._OncoSimul_test.hcbn)
## Done with the testing

## FIXME add similar testing for DiProg
## FIXME: change typeFitness to "modelType"



