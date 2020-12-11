Many calls return warnings and that is expected. In many calls I do not
care about the warnings and thus I do not test for them (I am testing
other things). However, there is no way to have testthat ignore them (see
https://github.com/r-lib/testthat/issues/772 and
https://stackoverflow.com/questions/41165707/force-testthat-to-ignore-warnings).
That forces me to wrap that in "supressWarning". This is work in progress.
