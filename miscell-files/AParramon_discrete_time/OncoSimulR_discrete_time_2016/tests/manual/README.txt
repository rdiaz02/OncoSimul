The tests under this directory are not run automatically when running "R
CMD check". They can be run by doing "test_dir()".

They are kept separate because they are long running, many are overkills,
and in a few it is possible some might fail even when there are no bugs
(since we are comparing expected outcomes from a distribution, and we
could end up in unlikely cases). In fact, you should expect some tests to
fail if you run them enough times (if they were to never fail, this could
mean that the tests are not actually sensitive enough).

Note that we will rarely want to use something like try_again (available
from the github version of testthat) or repeated tries (e.g.,
http://stackoverflow.com/questions/20770497/how-to-retry-a-statement-on-error)
as that could mean the condition is only true, say, 1/3 or 1/2 of the time
and that is NOT a proper test for these cases. However, we might want to
decrease the size of populations and allow failure in, say, less than 1 in
10 or similar. This way, we might run somewhat faster tests. 
