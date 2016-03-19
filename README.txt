cristinel.ababei@ndsu.edu
December 2008, Fargo ND


Synopsis
========
This is the REMNOC project: application REMapping for NOCs for fault 
tolerance. This is the code version that I used for experiments
submitted to RAW 2009. It works with only one application.

If you use this code or parts of it and would like to include a reference
to it, please use the following citation:
[1] C. Ababei and R. Katti, "Achieving Network-on-Chip fault tolerance by 
    adaptive remapping," IEEE Int. Parallel and Distributed Processing 
    Symposium (IPDPS), Reconfigurable Architectures Workshop (RAW), Rome, 
    Italy, May 2009.


Installation
============
First, edit the Makefile to reflect the location where you want
to compile and link remnoc. Then, just type:
> make


Remapping for fault tolerance
=============================
The new heuristic algorithm is described in the RAW 2009 paper.

The executable "remnoc" was created on RHEL (RedHat Enterprise Linux); 
in fact on Fedora 8. Current command options are:
  [-seed Int]
    Seed for the internal random number generator.
  [-failures_count Int]
    Number of failures to inject into the initial mapping.
  [-use_gui]
    If present, the gui will be launched too.
  [-remapping_algo heuristic | annealing]
   What algorithm to run: heuristic is the new algo, annealing is
   a in-house Simulated Annealing implementation for comparison purposes.

The gui is self-explanatory. If you want to save a .ps file with 
the currently displayed image, use "PstScript" button. The .ps file
will be created inside the "results" directory.


Examples
========
remnoc tests/vopd.app tests/vopd2.map -failures_count 1 -seed 3 -use_gui -remapping_algo heuristic
remnoc tests/vopd.app tests/vopd2.map -failures_count 1 -seed 3 -use_gui -remapping_algo annealing
remnoc tests/mpeg4.app tests/mpeg41.map -failures_count 1 -seed 3 -use_gui
remnoc tests/mpeg4.app tests/mpeg42.map -failures_count 1 -seed 3 -use_gui
remnoc tests/dsp_filter.app tests/dsp_filter1.map -failures_count 1 -seed 3 -use_gui
remnoc tests/dsp_filter.app tests/dsp_filter2.map -failures_count 1 -seed 3 -use_gui
remnoc tests/use_case2.app tests/use_case21.map -failures_count 1 -seed 3 -use_gui
remnoc tests/use_case2.app tests/use_case22.map -failures_count 1 -seed 3 -use_gui


Benchmarks
==========
See them in tests/.

Every benchmark is a set of at least two files. The first 
file, .app, represents the actual application graph. The second file,
and so on, .map, represents the initial mapping on a NOC of a 
certain size. I proposed these two formats because of the lack 
of any standards.

Example of .app file:
".application 2  <-- Application ID
.cores 6         <-- How many cores the application has
0 Mem            <-- Core ID and Core name. A new line for each core
1 ARM            <-- ...
2 Filter
3 IFFT
4 FFT
5 Display
.arcs 8          <-- How many arcs the application has
0 1 200          <-- Src Core ID, Des Core ID, Communication volume
1 2 600          <-- ...
2 1 600
1 5 200
2 3 200
3 2 200
2 4 200
4 2 200"

Example of .map file, used in conjunction with the .app above:
".mapping 2  <-- Mapping ID. We can have multiple .map files for an .app file
3 3          <-- N x N, the size of the NOC
0 0          <-- Core ID 0 is mapped to Tile ID 0
1 3          <-- Core ID 1 is mapped to Tile ID 3
2 6          <-- ...
3 7
4 4
5 1"

Tile IDs are assigned to tiles by counting all Tiles from left->right,
bottom->up. For example, the tile IDs of a 3x3 NOC architecture are:
 -----------
| 6 | 7 | 8 |
|-----------|
| 3 | 4 | 5 |
|-----------|
| 0 | 1 | 2 |
 -----------


Credits
=======
-- Vaughn Betz (while at Univ. of Toronto) developed much of the GUI;
-- Knuth's Stanford Graphbase: Hungarian Algorithm;


Copyright
=========
Copyright 2008 by Cristinel Ababei, cristinel.ababei@ndsu.edu
This Copyright notice applies to all files, called hereafter 
"The Software".
Permission to use, copy, and modify this software and its 
documentation is hereby granted only under the following 
terms and conditions.  Both the above copyright notice and 
this permission notice must appear in all copies of the 
software, derivative works or modified versions, and any 
portions thereof, and both notices must appear in supporting 
documentation.  Permission is granted only for non-commercial 
use.  For commercial use, please contact the authors.
This software may be distributed (but not offered for sale 
or transferred for compensation) to third parties, provided 
such third parties agree to abide by the terms and conditions
of this notice.
The Software is provided "as is", and the authors, the 
North Dakota State University (NDSU), as well as any and
all previous authors (of portions or modified portions of
the software) disclaim all warranties with regard to this 
software, including all implied warranties of merchantability
and fitness.  In no event shall the authors or NDSU or any and
all previous authors be liable for any special, direct, 
indirect, or consequential damages or any damages whatsoever
resulting from loss of use, data or profits, whether in an
action of contract, negligence or other tortious action,
arising out of or in connection with the use or performance
of this software.
