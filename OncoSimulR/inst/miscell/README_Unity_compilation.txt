What:
=====

This branch contains the C++ code of oncosimul using what is sometimes called "Unity builds":


- https://stackoverflow.com/a/33683596/3670562
- https://stackoverflow.com/a/3758580/3670562
- https://buffered.io/posts/the-magic-of-unity-builds/


In our case, the file unity_osimul.cpp includes the cpp files in the subdirectory UnityOncoSimul.


Pros:
=====

- Much faster install of the package because of much faster compilation. In Linux we go from about 5' for the install to about 2'. In Windows from about 47 to 57 minutes to about 17 to 18 (using, for optimization, -Os in both cases; with -O2 it takes about 19 to 20 minutes).

  Long installs were causing this to fail in the BioConductor servers (as times were often longer than 40 minutes)  



Cons:
====

See linked stackoverflow comments. Basically, changes in any one cpp file will lead to recompiles of the whole thing. Using ccache will not help much. But then, in Linux is about 2'.




Other: splitting file intervention.cpp
======================================

- File intervention.cpp has been split into several (intervention, intervention_2, intervention_private, ..., intervention_private_3) to make it possible, when not using the unity build, to compile it with optimization set to -O0 (to try to decrease speed). This is no longer an issue with the unity build, but might help in the future if we no longer use the unity build.


Further decreases in compilation time
=====================================

- Split the single cpp into two or four, and use MAKEFLAGS += -j2. Probably not worth it now.
- Note that decreases in compilation are not the only thing that matters, if the optimizations we leave out could make execution much faster.


Times (on my machine)
=====================

Using Windows 10 and Windows Server 2022 as virtual machines
(windows server: 1 CPU, windows 10 uses 2)

|           | Windows server | Windows 10 |
|-----------+----------------+------------|
| Non-unity |                |            |
| --------- |                |            |
|           |                |            |
| -Os       | 47'36"         |            |
|           |                |            |
| -O1       | 57'54"         |            |
|           |                |            |
| -O2       | 58'06"         | 57'04"     |
|           |                |            |
|           |                |            |
| Unity     |                |            |
| -----     |                |            |
|           |                |            |
| -O0       |                |            |
|           |                |            |
|           |                |            |
| -O1       |                | 18'35"     |
|           |                |            |
|           |                |            |
| -Os       | 18'00"         | 17'19"     |
|           | 16'51"         | 16'19"     |
|           |                |            |
| -O2       | 20'28"         | 19'42"     |
|           | 21'08"         | 20'46"     |
|           | 19'05"         |            |
|           |                |            |
| -O3       |                |            |
|           |                |            |
|           |                |            |




