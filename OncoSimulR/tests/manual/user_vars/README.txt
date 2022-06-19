The tests under this directory are the ones that check the user variable's functionality,
they have to check the whole x$other$userVarValues table, wich is evry big, in order to 
ensure that the variables get modified correctly, therefore, they are very long tests to execute.

A possible solution for this problem would be to modify executeRules in user_var.cpp so that it returns
true if a variable has been modified and false otherwise, and use this information in lines 1240-1253 of
BNB_nr.cpp to save the "currentTime" of the actualizations in a list that is included in OncoSimulIndiv's
return structure, when this is achieved the tests could then check the specific lines from x$other$userVarValues
where the times are the ones saved in this list, rather than checking the whole list, therefore the execution
time would be much lower.

Once this is made this tests should be moved to test.user_var.R in the testthat directory.