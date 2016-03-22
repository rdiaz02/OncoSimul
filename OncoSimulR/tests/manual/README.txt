The tests under this directory are not run automatically when running "R
CMD check". They can be run by doing "test_dir()".

They are kept separate because they are long running, many are overkills,
and in a few it is possible some might fail even when there are no bugs
(since we are comparing expected outcomes from a distribution, and we
could end up in unlikely cases).
