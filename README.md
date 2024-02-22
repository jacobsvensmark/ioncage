# ION and Charged Aerosol Growth Enhancement (ION-CAGE) model

The ION-Charged Aerosol Growth Enhancement (ION-CAGE) v0.9 code
constitutes a time-integration 0-dimensional model that evolves an
aerosol distribution through time while taking charge into account. The
model takes its initial and environmental conditions from an input-file,
and outputs snapshots with specified time intervals. If you use this
code, please cite [Svensmark et al. 
(2020)](https://doi.org/10.1029/2020EA001142).

## Installing, compiling and running the model

To install the model, simply clone this repository:

    > git clone https://github.com/jacobsvensmark/ioncage

In order to run the ION-CAGE model, the source code inside the `f90`
subdirectory should first be compiled using the `Makefile`:

    > make

This produces the executable `ioncage.out`. Note that this has only been
tested using the g95 FORTRAN compiler from OS X. If other compilers are
desired, the `Makefile` should be adjusted accordingly. If successful
the program may be run by typing the following

    > ./ioncage.out > outputfilename.dat

where `outputfilename.dat` is the data output file of the run (the
compiler may produce some warnings that can be ignored). Note that the
run is based on the input given in the `controlfile.txt` at runtime.

# Input parameters and inputfile

Most of the input for the model can be controlled by the
`controlfile.txt` file, and it must be placed in the same directory as
the compiled executable file produced by the `Makefile`. The
`controlfile.txt` is split up into a number of paragraphs:

-   In the **Total number of particle nodes** top paragraph of the
    `controlfile.txt`, the number of volume nodes can be set, i.e. the
    number of log-spaced aerosol sizes.

-   The second paragraph of `controlfile.txt` it is possible to switch
    on (`FLAG=1`) and off (`FLAG=0`)
    different modules i.e. different terms of the general dynamic
    equation that the code runs. Each module corresponds to
    contributions to $dN/dt$ via different mechanisms (loss, coagulation
    etc), and the code for each module can be found in the
    `dNdtmodule.f90` file.

-   The **Continue from previous run** paragraph is concerned with the
    possibility of taking an end-state of a previous run as initial
    conditions for the current run. The name of the previous run
    data-file should be provided, and the FLAG set to 1 to actually use
    it.

-   In the **Parameters** paragraph several parameters may be set for
    the course of the run. First the nucleation rate of small aerosols
    in neutral or charged state. Note that the nucleation of of a
    charged aerosol removes an ion of the same charge. Then the
    production rate for monomers (neutral and charged), the ion-ion
    recombination coefficient, neutral monomer density and ion
    densities, loss terms, masses (should match the used interaction
    coefficients).

-   In the **Initial concentrations** paragraph initial ion and sulfuric
    acid concentrations are set. The `fixinitials` keyword has three
    flags separated by comma, controlling whether the initial
    concentrations of the positive ions, negative ions and sulfuric acid
    are kept constant throughout the entire run (in that order). Initial
    values are overwritten if `Loadfileflag=1`, however if the fix
    `fixinitials` are set to 1, the choice of fixed ions and sulfuric
    acid are not overwritten. Also the initial values may be kept
    constant by setting the `fixinitials=1`.

-   The **Time related parameters** control start time (usually 0),
    `End time` i.e. length of run, step length $\Delta t$ for the
    integration and the time between each snapshot printed to the output
    file.

-   The **Size related parameters** paragraph contains parameters
    controlling the size range of the nodes in the simulation, as well
    as critical diameters for neutral and charged nucleated aerosols.

-   Finally, the **Interaction Coefficients** paragraph contains a list
    of paths to datafiles with interaction coefficients that are loaded
    for each type of interaction used in the GDE. (See subroutines
    `readcoefftwocolumn` and `readcoeffmatrix` in file
    `ioioncagemodule.f90` for an explanation on the format assumed for
    the coefficient data files). For each file, a number of header lines
    to skip before reading data can be specified. **NOTE:** It is
    important to be aware exactly how these are used in the `main.f90`
    and `dNdtmodule.f90`, as many different approaches are possible
    here.

Each parameter has a dedicated keyword followed by a value, which should
be self explanatory. Lines that do not start with one of these keywords
are ignored. Note that the executable does **not** need to be recompiled
to register changes in `controlfile.txt`.

## Units

Unless otherwise stated, all calculations internally and all output are
in standard SI units. Note that several input parameters in the
`controlfile.txt` are not necessarily SI units, but converted to SI by
the `readcontrolfile` subroutine.

## CPU-Time vs. Precision

The model run speed varies significantly depending on the calculation.
The number of nodes and temporal resolution can be adjusted to decrease
computation time at the cost of precision inside the `controlfile.txt`.
Especially the coagulation subroutine is sensitive to the number of
nodes as it is an $N\times N$ process. The coagulation subroutine may be
disabled by setting the corresponding FLAG inside `controlfile.f90`.

## Visualizing output

A couple of python scripts inside of the `py` subdirectory have been
included for reading and visualizing model output.

## File overview

Files for running the model, all coded in FORTRAN include:

-   `main.f90`: The top level program i.e. the main program for the
    model. It is here the Runge Kutta time integration loop of the
    General Dynamics Equation is found.

-   `dNdtmodule.f90`: Module containing subroutines for the calculation
    of coagulation, condensation, ion-induced condensation, nucleation,
    charge exchange and loss terms as controlled in the
    `controlfile.txt`. Called from main.f90.

-   `ioioncagemodule.f90`: Subroutines related to I/O for the model.
    `readcoefmatrix` and `readcoeftwocolumn` are used to read
    interaction coefficients from a file (see instructions in the source
    code to use your own tables). `readcontrolfile` reads all parameters
    set in controlfile and converts to SI units.
    `printcontrolfiletooutput` is used to print the contents of the used
    controlfile into the output. `printsnapshot` prints a snapshot of
    the simulation to the output. Finally `loadstate` loads the latest
    snapshot of an output data from a previous run, and continues
    simulation from that state if corresponding FLAG is set in
    `controlfile.txt`.

-   `mathmodule.f90`: This modile contains two interpolation subroutines
    `linearinterpolation` and `cubicinterpolation2d` used in the process
    of loading interaction coefficients from the `readcoeftwocolumn` and
    `readcoefmatrix` subroutines. Both will return an error if grid is
    out of the interpolation range. Therefore, remember to check if
    diameter range of the model is compatible with (contained in) the
    range of diameters provided in these files.

-   `iogeneral.f90`: Contains `filelines` which counts the number of
    lines in a text file, and `findnl` which locates newline character
    in a string.

-   `precisiontype.f90`: FORTRAN module specifying precision of floats.

-   `controlfile.txt` Inputfile controlling parameters and other types
    of setup for the model at runtime. See below for in-depth
    explanation of the controlfile contents.

Also included are a number of interaction coefficient tables for
condensation and coagulation. Coagulation coefficient tables are
provided for aerosols of density 1200kg/m3. Condensation
coefficients are provided for four different masses of monomers, namely
100AMU, 150AMU, 225AMU and 350AMU, which can be selected
in the controlfile.
