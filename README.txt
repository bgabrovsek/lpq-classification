Classification of knots in lens spaces
Created by Bostjan Gabrovsek 31/1/13, bostjan.gabrovsek@fmf.uni-lj.si.

FILES:

main.cpp               - main function (starting execution),
knot.h                 - knot class and Gauss code manipulation functions,
link.h                 - link class (used by HOMLFY.h), link manipulation functions,
KBSM.h                 - KBSM of knots in the solid torus,
HOMFLY.h               - HOMFLY skedin module of knots in the solid torus,
lens.h                 - manipulation of links in L(p,q): SL move, winding number, HSM(L(p,q)), KBSM(L(p,q)),
multivariate_laurent.h - multivariable laurent polyomial library for up to 256 variables,
global_vars.h          - global variables and definitions used throughout the project space (file names, constants)
numbers.h              - classes of 128-bit, 256-bit, and 512-bit integers,
in_out.h               - input, output, and print functions
algorithms.h           - main algorithms: partitioning knots by invariants, reducing knots, classification algorithms,...
classification_data.h  - classification class that includes the information if knots are prime or duplicates,
timing.h               - timer class (for measuring speed of execution),
tree.h                 - tree class, used for manipulating arrow diagrams,
common.h               - various constants and macros (string parsing, bit-wise operations, memory functions, compare functions),
randomized.h           - random Reidemeister walk,
knot_preprocessor_directives.h - preprocessor directives (sets architecture depending on the maximal number of crossings).

results.zip            - output copies (i.e. results) of the classification algorithm


COMPILATION

For optimal performance, run each step with the program compiled with the suggested maximal number of crossings.
Running the code with a higher number of crossings than necessary will take more time.

INPUT ARGUMENTS

Only one argument can be used at a time.

-g
 
Generates all ralizable Gauss codes up to N crossings, saves the GC to file,
calculates the HOMLFY skein modules partitions,
sets first knot in each partition as prime
stores the classification to file
(for this option code should be compiled with MAX >= 12).
 
-c
 
Tries to reduce each non-prime knot using R-moves, flypes, and meridional rotations
(for this option code should be compiled with MAX >= 12).
 
-h
Calculates HOMLFY skein module generators in terms of solid torus HOMLFY skein module generators
(for this option code should be compiled with MAX >= 220).
May need a lot of time to execute, the step can be avoided by copying lpq-homfly*.txt to the working directory.
 
-i
Initialized classification of knots in L(p,q)
(for this option code should be compiled with MAX >= 12).
 
-l
Tries to reduce non-prime knots in L(p,q) using R-moves, flypes, and meridional rotations
(for this option code should be compiled with MAX >= 127).
 
-m
Manualy reduce L(p,q) knots by R-moves, flypes, meridional rotations, and slides
(for this option code should be compiled with MAX => 63).
 
-pt
Prints the classification table(s)
(for this option code should be compiled with MAX >= 12).

-pf
Prints knots that share HSM or KBSM
(for this option code should be compiled with MAX >= 12).
 
-ph
Prints HSM of knots
(for this option code should be compiled with MAX >= 12).

-pm
Prints KBSM of knots
(for this option code should be compiled with MAX >= 12).

-r
Random walks of unclassified knots
(for this option code should be compiled with MAX => 31).




