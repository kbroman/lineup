## Version 0.44, 2024-07-15

- Fixed link in documentation


## Version 0.43-1, 2023-02-19

- Revise CITATION file to use `bibentry()` rather than `citEntry()`


## Version 0.42, 2022-07-10

- Include <limits.h> in a C file, for an upcoming change in R


## Version 0.40, 2021-11-03

- Fixed lineupversion() to handle a case like "0.40".


## Version 0.38-3, 2020-11-16

- Fixed use of class(), for example using inherits(object, "blah")
  rather than "blah" %in% class(object).

- Fixed potential problem in documentation, since plot() has moved
  from the graphics package to base.

- Converted documentation to use markdown


## Version 0.37-11, 2019-03-06

- Fixed a bug in distee()

- Fixed test for corbetw2mat(), due to upcoming change in
  sample(), in R 3.6


## Version 0.37-10, 2018-03-23

- Add ByteCompile to the DESCRIPTION file


## Version 0.37-9, 2017-05-19

- Small changes to avoid Note about "R_registerRoutines"


## Version 0.37, 2015-07-20

- Add a vignette and example data.

- calc.locallod can do calculations in parallel


## Version 0.36, 2014-05-02

- Convert documentation to use roxygen2


## Version 0.35, 2012-10-25

- Only minor changes


## Version 0.34, 2012-10-15

- In plot2dist, make default smoothScatter=FALSE


## Version 0.33, 2012-08-14

- Added function combinedist, for combining distance matrices


## Version 0.32, 2012-08-08

- Added addcovar and intcovar arguments to calc.locallod.


## Version 0.31, 2012-03-02

- Added NAMESPACE.

- Revised lineupversion() to use packageVersion().


## Version 0.30, 2011-05-11

- This is a new package to contain tools for lining up two data
  matrices, to detect and correct sample mix-ups.  The key examples
  are gene expression microarray data on a large set of samples on
  two tissues, or between genotype and gene expression.
