The source code and makefile are in this directory, the test inputs are in the
ex subdirectory, and the results are in results.txt in the Table 1 order.

The tests are built and run as follows:
% make test1
% test1 ex/cube.vtk ex/glacier.vtk ex/cube.vtk
and likewise for tests 2, 3, 4, 6, and 8.  Each test outputs a one line result
(except that test 6 outputs multiple lines) followed by statistics in the Table
2 format.

Tests 5 and 7 are combined and have a different format
% make test57a
% test57a ex/frustum ex/tworooms
and
% make test57b
% test57b ex/plus ex/lattice-room
The output is verbose.  Test 5 concludes with a printout of its statistics,
test 7 begins, and finally test 7 ends with a printout of its statistics.

Tests 1-4 and 8 can be run with arbitrary inputs, test 6 has no input, and
tests57 is restricted to the above two inputs.

To replace PDD with exact evaluation of rational predicates, remove the # from
#-DEXACTRAT in line 1 of makefile, call % make clean, and make the test.

The default version of test 6 uses the perturbation identity detection algorithm
and bivariate predicates (n = 2 in Table 4).  The univariate (n = 1) test is
obtained by uncommenting line 25 in angle.h:
//#define HOMOTOPY_IDENTITY
The residue algorithm versions with n = 1 and n = 2 are obtained by uncommenting
line 23
//#define RESIDUE_IDENTITY1
and line 24
//#define RESIDUE_IDENTITY2
respectively.

To run the table 4 variants of test 6, copy acp-test6.h to acp.h, % make clean,
% make test6. After the test, copy acp-standard.h to acp.h


