# TGMG: A templated geometric multigrid library in C++
TGMG is a software library for solving a class of [elliptic partial differential
equations](https://en.wikipedia.org/wiki/Elliptic_partial_differential_equation)
on two-dimensional rectangular grids using the [multigrid
method](https://en.wikipedia.org/wiki/Multigrid_method). A defining
characteristic of the library is that it makes use of C++ templates.
This allows it to work on a range of data types, and also allows
performance-critical parts of the solution method to be inlined during
compilation.

The library makes use of [OpenMP](https://www.openmp.org) for multithreading,
and self-tunes the number of threads used in different parts of the algorithm.
It works on any m &times; n rectangular grid, using either periodic or
non-periodic boundary conditions.

The library will only solve linear systems where the non-zero terms connect
each grid point to its 3 \times; 3 square of orthogonal and diagonal neighbors,
but this includes a wide-range of systems arising in finite-difference and
finite-element discretizations.

The library is [header-only](https://en.wikipedia.org/wiki/Header-only). To use
the library, the user must write a problem class that encapsulates all of the
details of the linear system to be solved. Several examples are provided. This
documentation assumes that the reader is familar with the principles behind
multigrid methods; for an introduction, see the textbooks by Briggs *et al.*
[1] and Demmel [2].

## Compiling the code examples
The code is written in C++ and has been tested on Linux, MacOS, and Windows
(via [Cygwin](https://www.cygwin.com)). It uses [Perl](http://www.perl.org)
during the compilation procedure. The following documentation assumes you are
familiar with the Linux/Mac/Cygwin
[command-line interface](https://en.wikipedia.org/wiki/Command-line_interface).

The code does not have any dependencies, although it can output data in a
binary format that can be read by the freeware plotting program
[Gnuplot](http://www.gnuplot.info).

To compile the code it is necessary to create a common configuration file
called **config.mk** in the parent directory where the **tgmg** repository is
stored. Several templates are provided in the **config** directory. To use,
copy one of the templates into the parent directory. From the **tgmg**
directory, on a Linux computer, type
```Shell
cp config/config.mk.linux ../config.mk
```
On a Mac using GCC 10 installed via [MacPorts](http://www.macports.org), type
```Shell
cp config/config.mk.mac_mp ../config.mk
```
On a Mac using GCC installed via [Homebrew](http://brew.sh), type
```Shell
cp config/config.mk.mac_hb ../config.mk
```
On a Windows computer with Cygwin installed, type
```Shell
cp config/config.mk.win_cw ../config.mk
```
After this, the code examples can be compiled by typing
```Shell
make
```

## Finite-difference example
The example program **pfd_test.cc** demonstrates using the library to solve the
Poisson equation &nabla;&sup2;u=f on a two-dimensional domain [-1,1]&sup2; with
zero Dirichlet boundary conditions on all sides. The solution is discretized on
a 1025 &times; 1025 grid, and the Poisson equation is discretized using the
standard [five-point centered-finite difference
stencil](https://en.wikipedia.org/wiki/Discrete_Laplace_operator) for the
Laplacian. The corresponding linear system Az=b is specified in the
**poisson_fd** class, which is described in more detail in the following
section.

The first part of the **pfd_test.cc** program creates a **poisson_fd** class
called **pfd** with the specified grid size and grid dimensions:
```C++
    // Grid dimensions
    const int m=1025,n=1025;

    // Physical dimensions of the rectangle to solve on
    const double ax=-1,bx=1,ay=ax,by=bx;
    poisson_fd pfd(m,n,ax,bx,ay,by);
```
After this the multigrid solver is created. It is a template that is
instantiated on three classes: (1) the problem class, (2) the data type for
entries in A, and (3) the data type for entries in b.
```C++
    // Initialize the multigrid solver, set the verbosity to maximum, and set
    // up the multigrid hierarchy
    tgmg<poisson_fd,double,double> mg(pfd,pfd.b,pfd.z);
```
TGMG can has four different levels of messages, which is controlled by the
**verbose** class member. They are: (0) silent operation, (1) error messages
only, (2) a summary of convergence information for the solution, (3) detailed
residual messages per V-cycle. In this example, the level is set to three. In
addition, a call to the **setup()** function is made, which calculates the
linear systems on the coarser grids:
```C++
    mg.verbose=3;
    mg.setup();
```
The source term is then initialized as a Gaussian, and the solution array is
set to zero, after which the linear system is solved using multigrid V-cycles:
```C++
    // Set up the solution and source arrays
    mg.clear_z();
    pfd.gaussian_source_term(0,0,1,1);

    // Solve using multigrid V-cycles
    mg.solve_v_cycle();
```
Finally, the solution and source terms are saved in a binary format that can be
read by Gnuplot:
```C++
    // Output the solutions in a format that can be read by Gnuplot using
    // either of the following two commands:
    // 
    // 1. For color map plot
    // plot [-1:1] [-1:1] 'filename' matrix binary with image
    //
    // 2. For 3D plot
    // splot [-1:1] [-1:1] 'filename' matrix binary
    const double dx=(bx-ax)/(m-1),dy=(by-ay)/(n-1);
    mg.output_b("src.0",ax,dx,ay,dy);
    mg.output_z("sol.0",ax,dx,ay,dy);
```
When the program is run using eight threads it produces the following output:
```
# Top grid level: (1025,1025) [odd,odd] {8,8,-,8}
# Grid level  0 : (513,513) [odd,odd] {8,8,8,8}
# Grid level  1 : (257,257) [odd,odd] {8,8,8,8}
# Grid level  2 : (129,129) [odd,odd] {8,8,1,1}
# Grid level  3 : (65,65) [odd,odd] {1,1,1,1}
# Grid level  4 : (33,33) [odd,odd] {1,1,1,1}
# Grid level  5 : (17,17) [odd,odd] {1,1,1,1}
# Grid level  6 : (9,9) [odd,odd] {1,1,1,1}
# Grid level  7 : (5,5) [odd,odd] {1,-,1,1}
Iteration 4, residual 5.18612e-10
Iteration 8, residual 1.95001e-18
8 iters, res 0.0468185->1.95001e-18, 2.04755 digits per iter
```
The first few lines describe the hierarchy of coarser grids that TGMG has
created. The grids are recursively constructed, with an m &times; n grid being
coarsened to &lfloor;(m+1)/2&rfloor; &times; &lfloor;(n+1)/2&rfloor; grid.
For non-periodic grids, using problem dimensions of 2<sup>k</sup>+1 gives a
convenient implementation [2], since every grid in the hierarchy has an odd
number of points. Restriction on an odd grid size works simply since coarsening
works by removing every other grid point (while retaining the end points).
However, TGMG can also coarsen grid dimensions with even sizes, allowing it to
work on any grid.

The numbers in the curly brackets correspond to the threads that TGMG has
chosen to perform the different operations on that grid. The four numbers
correspond to (1) Gauss–Seidel sweeps, (2) restriction, (3) interpolation, and
(4) clearing arrays. Once the grid becomes too small, it is no longer
advantageous to use threads, and the library switches to single-threaded
operation.

The library then prints out the square residual per grid point as the
iterations progress. After eight iterations the error threshold has been
reached, and the library prints a message about the overall reduction in
the residual. The code estimates the number of digits of accuracy gained per
iteration, but this is only approximate since once the algorithm converges, the
residual may be dominated by rounding error.

The images below show the source term (left) and the solution (right).

![Plots of Gaussian source term (left) and solution to the Poisson equation](https://people.seas.harvard.edu/~chr/tgmg/pfdt.png)

## Structure of the problem class
The linear system in the previous example is specified in the **poisson_fd**
class, contained in the **poisson_fd.cc** and **poisson_fd.hh** source code
files. This class contains the following class data members required by TGMG:

- **m** and **n**, the grid dimensions,
- **x_prd** and **y_prd**, boolean values setting the periodicity in the two
  directions,
- **gs_mode**, the type of Gauss–Seidel sweep to use in the solver,
- **acc**, the solution tolerance on square error per grid point,
- **z**, the solution vector array,
- **b**, the source vector array.

In addition, the class specifies several functions that describe the linear
system to be solved. The functions are of the horizontal index **i** and the
overall grid index **ij** (equal to **i+j&times;m**). The vertical index
is not given, because it is possible to infer it using simple calculations from
the other two.

Nine functions set the terms in the linear system:

|**a_ul(i,ij)**|**a_uc(i,ij)**|**a_ur(i,ij)**|
|--------------|--------------|--------------|
|**a_cl(i,ij)**|**a_cc(i,ij)**|**a_cr(i,ij)**|
|**a_dl(i,ij)**|**a_dc(i,ij)**|**a_dr(i,ij)**|

For a given grid point, these functions return the matrix entries corresponding
to the 3 &times; 3 set of neighboring grid points.

Two further redundant functions are provided:

- **inv_cc(i,ij,v)** returns the division of the value v by the diagonal entry
  of the linear system at grid point (i,j),
- **mul_a(i,ij)** returns the result of the matrix multiplication A'x at grid
  point (i,j), where A' contains all off-diagonal entries in A.

Both of these functions are frequently called in the Gauss–Seidel sweeps, and
they are separately specified to increase computational efficiency. The
**poisson_fd** class contains other members that are specific to the
finite-difference implementation (*e.g.* hcc, hxx). In general, these problem
classes can contain many other types of functions and data.

## Other test problems
The program **pfd_time** measures the time for the TGMG library to solve the
**poisson_fd** problem 500 times, using Gaussian source terms that slowly
change for each test case. On a 1025 &times; 1025 grid, using eight threads on
an iMac with 3.6 GHz 8-core Intel Core i9, the code takes 43 ms per solution,
with an average of six V-cycles per solve.

The program **pfem_test** demonstrates solving a Poisson equation using
a finite-element discretization based on bilinear equations. It makes use of
the **poisson_fem** problem class described in **poisson_fem.hh** and
**poisson_fem.cc**.

The program **vpfem_test** solves the equation -&nabla;&middot;(c&nabla;u)=f
where c is a spatially-varying field. Equations like this occur in many
physical problems, such as when simulating fluids with varying density [3,4,5],
or simulating porous media flow with varying permeability [6].

If c has large variations in size over the domain, then V-cycles become
inefficient. The program **mpcg_test** uses TGMG as a component in the
multigrid preconditioned conjugate gradient algorithm [7], which gives better
convergence.

## Applications of the library
This library has been developed by Chris Rycroft with feedback from the
Rycroft group and collaborators, and it has been used in a number of scientific
projects and publications.

It has been used in a sequence of papers that study bulk metallic glasses
[8,9,10], a new type of alloy under consideration for a variety of technological
applications. It has been used in the incompressible reference map technique [5],
a new numerical method for fluid–structure interaction. In has also been use
to inforce incompressibility constraints in models of porous media flow [6,11].

## Contact
For questions about the code, contact [Chris Rycroft](http://seas.harvard.edu/~chr/).

## Acknowledgements
This work has been partially supported by the National Science Foundation under
Grant Nos. DMR-1409560 and DMS-1753203, and by the Applied Mathematics Program
of the U.S. DOE Office of Science Advanced Scientific Computing Research under
contract number DE-AC02-05CH11231.

## Bibliography
1. James W. Demmel, *Applied Numerical Linear Algebra*, SIAM (1997).
   [doi:10.1137/1.9781611971446](https://doi.org/10.1137/1.9781611971446)

2. William L. Briggs, Van Emden Henson, and Steve F. McCormick, *A Multigrid
   Tutorial, Second Edition*, SIAM (2000).
   [doi:10.1137/1.9780898719505](https://doi.org/10.1137/1.9780898719505)

3. Mark Sussman, Ann S. Almgren, John B. Bell, Phillip Colella, Louis H.
   Howell, and Michael L. Welcome, *An adaptive level set approach for
   incompressible two-phase flows*, J. Comput. Phys. **148**, 81–124 (1999).
   [doi:10.1006/jcph.1998.6106](https://doi.org/10.1006/jcph.1998.6106)

4. Jiun-Der Yu, Shinri Sakai, and James A. Sethian, *A coupled level set
   projection method applied to ink jet simulation*, Interface Free Bound.
   **5**, 459–482 (2003).
   [doi:10.4171/IFB/87](https://doi.org/10.4171/IFB/87)

5. Chris H. Rycroft, Chen-Hung Wu, Yue Yu, and Ken Kamrin, *Reference map
   technique for incompressible fluid–structure interaction*, J. Fluid Mech.
   **898**, A9 (2020).
   [doi:10.1017/jfm.2020.353](https://doi.org/10.1017/jfm.2020.353)

6. Nicholas J. Derr, David C. Fronk, Christoph A. Weber, Amala Mahadevan, Chris
   H. Rycroft, and L. Mahadevan, *Flow-driven branching in a frangible porous
   medium*, [arXiv:2007.02997](https://arxiv.org/abs/2007.02997) (2020).

7. O. Tatebe, *The multigrid preconditioned conjugate gradient method*, in 6th
   Copper Mountain Conference on Multigrid Methods, Copper Mountain, CO, April
   4–9, 1993. [NASA Technical Reports Server 19940017009](https://ntrs.nasa.gov/search.jsp?R=19940017009)

8. Chris H. Rycroft, Yi Sui, and Eran Bouchbinder, *An Eulerian projection
   method for quasi-static elastoplasticity*, J. Comput. Phys. **300**,
   136–166 (2015).
   [doi:10.1016/j.jcp.2015.06.046](https://doi.org/10.1016/j.jcp.2015.06.046)

9. Manish Vasoya, Chris H. Rycroft, and Eran Bouchbinder, *Notch fracture
   toughness of glasses: Rate, age and geometry dependence*, Phys. Rev. Applied
   **6**, 024008 (2016).
   [doi:10.1103/PhysRevApplied.6.024008](https://doi.org/10.1103/PhysRevApplied.6.024008)

10. Adam R. Hinkle, Chris H. Rycroft, Michael D. Shields, and Michael L. Falk,
   *Coarse graining atomistic simulations of plastically deforming amorphous
   solids*, Phys. Rev. E. **95**, 053001 (2017).
   [doi:10.1103/PhysRevE.95.053001](https://doi.org/10.1103/PhysRevE.95.053001)

11. Christoph A. Weber, Chris H. Rycroft, and L. Mahadevan, *Differential
   activity-driven instabilities in biphasic active matter*, Phys. Rev. Lett.
   **120**, 248003 (2018).
   [doi:10.1103/PhysRevLett.120.248003](https://doi.org/10.1103/PhysRevLett.120.248003)
