.. _case:

Case files
==========

NekRS simulation requires a number of case definition files which are described in this page.
An overview of these are presented in the image below .

.. _fig:case_overview:

.. figure:: ../_static/img/overview.svg
   :align: center
   :figclass: align-center
   :alt: An overview of nekRS case files

There are a minimum of three files required to run a case:

* **Parameter file**, with ``.par`` extension. This sets simulation parameters used by the case and can be modified between runs.
* **User-defined host file**, with ``.udf`` extension. This is used to set specific equations of the case, initial/boundary conditions, data outputs and other user definable behaviour.
* **Mesh file**, with ``.re2`` extension. This file defines the geometry of the case.

Some optional files can also be included:

* **User-defined okl file**, with ``.oudf`` extension. This file conventionally has all the user defined device kernels. It must be included in :ref:`okl_block` of ``.udf`` file.
* **Legacy Nek5000 user file**, with ``.usr`` extension. This file allows usage of *Nek5000*, legacy, Fortran 77 user routines.

The case name is the default prefix applied to these files - for instance, a complete input description with a case name of "eddy" would be given by the files ``eddy.par``, ``eddy.re2``, ``eddy.udf``.
Optionally, the user can also define the names of ``.udf``, ``.oudf`` and ``.usr`` in the :ref:`sec:generalpars` section and name of ``.re2`` file in :ref:`sec:meshpars` section of ``.par`` file. 

The following sections describe the structure and syntax for each of these files for a general case.


.. _parameter_file:

Parameter File (.par)
---------------------

.. tip::

   The user can access the manual containing specifications of the parameter file using ``nrsman`` as follows

   .. code-block:: bash

      nrsman par

Most information about the problem setup is defined in the parameter (``.par``) file.
This file is organized in a number of **sections**, each with a number of **keys**.
Values are assigned to these keys in order to control the simulation settings.

The general structure of the ``.par`` file is as follows, where ``FOO`` and ``BAR`` are both section names, with a number of (key, value) pairs.

.. code-block:: ini

  [FOO]
    key = value
    baz = bat

  [BAR]
    alpha = beta
    gamma = delta + keyword=value + ... 

The valid sections for setting up a *NekRS* simulation are:

* ``GENERAL``: generic settings for the simulation (mandatory, see :ref:`General Parameters<sec:generalpars>`)
* ``OCCA``: backend :term:`OCCA` device settings (optional, see :ref:`OCCA section<sec:occa>`)
* ``PROBLEMTYPE``: settings for the governing equations (optional, see :ref:`Problem Type<sec:problemtype>`)
* ``MESH``: mesh related settings (optional, see :ref:`Mesh Parameters<sec:meshpars>`)
* ``GEOM``: **TODO**
* Field Settings (see :ref:`Field Settings<sec:field_settings>`)

  * ``FLUID VELOCITY``: settings for the velocity solver

  * ``FLUID PRESSURE``: settings for the pressure solver 

  * ``SCALAR``: default scalar settings

    * ``SCALAR FOO``: settings for the ``FOO`` scalar

* ``BOOMERAMG``: settings for the Hypre's :term:`AMG` solver
* ``NEKNEK``: settings for the *NekNek* module in *NekRS* (see :ref:`NekNek Parameters <sec:neknekpars>`)
* ``CVODE``: settings for the CVODE solver (see :ref:`CVODE Parameters <sec:cvodepars>`)
  
.. note::

  - Section name and key/value pairs are treated as case insensitive
  - Values enclosed within quotes maintain case sensitivity
  - Values prefixed with 'env::' are interpreted as references to environment variables

.. _sec:user_section:

User Sections
""""""""""""""""""

The user also has the option to specify additional sections to define custom control keys in ``.par`` file.
These sections must be declared at the top of the ``.par`` file using ``userSections`` key as shown in the below example


.. code-block:: ini

   userSections = CASEDATA

   ...

   [CASEDATA]
   key = value


.. _sec:generalpars:

General Parameters
""""""""""""""""""

.. _tab:generalparams:

.. csv-table:: ``GENERAL`` keys in the ``.par`` file
   :widths: 20,20,60
   :header: Key, Value(s), Description/Note(s)/Default Value

   ``polynomialOrder``,``<int>``, "``polynomialOrder`` > 10 is currently not supported"
   ``dealiasing``,``true`` / ``false``, "Enables/disables over-integration of convective term |br| Default = ``true``"
   ``cubaturePolynomialOrder``,``<int>``, "Polynomial order of ``dealiasing`` |br| Default = 3/2*(``polynomialOrder`` +1)-1"
   ``verbose``,``true`` / ``false``, "``true`` instructs *NekRS* to print detailed diagnostics to *logfile* |br| Default = ``false``"
   ``redirectOutputTo``,``<string>``,"String entry for the name of the *logfile* to direct *NekRS* output"
   ``startFrom``,"``<string>`` |br| ``+ time=<float>`` |br| ``+ x`` |br| ``+ u`` |br| ``+ s or s00 s01 s02 ...`` |br| ``+ int``", "Restart from specified ``<string>`` file |br| reset ``time`` to specified value |br| read mesh coordinates |br| read velocity |br| read all scalar or specified scalars |br| interpolate solution (useful if mesh coordinates are different)" 
   ``timeStepper``,``tombo1`` / ``tombo2`` / ``tombo3``," Order of time discretization for BDFk/EXTk scheme |br| Default = ``tombo2``"
   ``stopAt``,``numSteps`` / ``endTime`` / ``elapsedTime``, "stop criterion |br| Default = ``numSteps``"
   ``numSteps``,``<int>``, "Number of simulation time steps"
   ``endTime``,``<float>``,"Simulation end time"
   ``elapsedTime``,``<float>``,"Simulation time in wall clock minutes"
   ``dt``,``<float>`` |br| ``+ targetCFL = <float>`` |br| ``+ max = <float>`` |br| ``+ initial = <float>`` , "Time step size |br| adjust ``dt`` to match ``targetCFL`` |br| max limit of ``dt`` |br| Initial ``dt`` "
   ``advectionSubCyclingSteps``,``<int>``,"Number of OIFS sub-steps for advection |br| Default = ``0`` (OIFS turned off)"
   ``constFlowRate``,"``meanVelocity = <float>`` |br| ``meanVolumetricFlow = <float>`` |br| ``+ direction = <X,Y,Z>``","Specifies constant flow velocity |br| Specifies constant volumetric flow rate |br| Specifies flow direction" 
   ``scalars``,"``<string>, <string> ...``","Name of scalar fields to be solved"
   ``checkPointEngine``,``<string>`` |br| ``nek`` / ``adios``,"Specifies engine to write field files |br| Default = ``nek``"
   ``checkPointPrecision``,``<int>`` |br| ``32`` / ``64``,"Specifies precision of field files |br| Default = ``32``"
   ``checkPointControl``,``steps`` / ``simulationTime``,"Specifies check point frequency control type |br| Default = ``steps``"
   ``checkPointInterval``,``<int>`` / ``<float>`` |br| 0 |br| -1, "Specifies check point frequency (``<int>`` for ``steps`` / ``<float>`` for ``simulationTime``) |br| ``0`` implies at end of simulation |br| ``-1`` disables checkpointing" 
   ``udf``,"``''<string>''``","Optional name of user-defined host function file |br| Default is ``<case>.udf``"
   ``oudf``,"``''<string>''``","Optional name of user-defined OCCA kernel function file |br| As a default *NekRS* expects these are defined in :ref:`OKL block <okl_block>` in ``.udf`` file"
   ``usr``,"``''<string>''``","Optional name of user-defined legacy *Nek5000* (fortran) function file |br| Default is ``<case>.usr``"
   ``regularization``,"","Specifies regularization options for all fields |br| See :ref:`common field settings<sec:common_settings>` for details"

.. _sec:occa:

OCCA Parameters
""""""""""""""""
.. _tab:occaparams:

.. csv-table:: ``OCCA`` keys in the ``.par`` file
   :widths: 20,20,60
   :header: Key, Value(s), Description/Note(s)/Default Value

   ``backend``, |br| ``SERIAL`` / |br| ``CUDA`` / |br| ``HIP`` /|br| ``DPCPP``,"Specifies the *device* for JIT compilation. Default is defined ``$NEKRS_HOME/nekrs.conf`` |br| CPU |br| NVIDIA GPU (CUDA) |br| AMD GPU (HIP) |br| Intel GPU (oneAPI)"
   ``deviceNumber``,``<int>`` |br| ``LOCAL-RANK``,"Default is ``LOCAL-RANK``"
   ``platformNumber``,``<int>``, "Only used by ``DPCPP`` |br| Default is ``0``"

.. _sec:problemtype:

Problem Type Parameters
""""""""""""""""""""""""""
.. _tab:problemparams:

.. csv-table:: ``PROBLEMTYPE`` keys in the ``.par`` file
   :widths: 20,20,60
   :header: Key, Value(s), Description/Note(s)/Default Value

   ``equation``,``stokes`` |br| ``navierStokes`` |br| ``+ variableViscosity``, "Stokes solver |br| Navier-Stokes solver |br| uses stress formulation (required for spatially varying viscosity)"

.. _sec:meshpars:

Mesh Parameters
""""""""""""""""
.. _tab:meshparams:

.. csv-table:: ``MESH`` keys in the ``.par`` file
   :widths: 20,20,60
   :header: Key, Value(s), Description/Note(s)/Default Value

   ``partitioner``,``rbc`` / ``rsb`` / ``rbc+rsb``,"Specifies mesh partitioner |br| Default = ``rbc+rsb`` "
   ``boundaryIDMap``,"``<int>, <int>, ...``", "Map mesh boundary ids to 1,2,3,... |br| See :ref:`boundary conditions<boundary_conditions>` for details"
   ``boundaryIDMapFluid``,"``<int>, <int>, ...``", "Required for conjugate heat transfer cases |br| See :ref:`boundary conditions<boundary_conditions>` for details"
   ``connectivityTol``,"``<float>``","Specifies mesh tolerance for partitioner |br| Default = ``0.2``"
   ``file``,"``''<string>''``","Optional name of mesh (``.re2``) file |br| Default is ``<case>.re2``"


.. _sec:field_settings:

Field Settings
"""""""""""""""""""""

The sections for specific fields, including velocity (``FLUID VELOCITY``), pressure (``FLUID PRESSURE``) and scalars (``SCALAR`` or ``SCALAR FOO``) contain keys to describe linear solver setting for the corresponding field.
Most of the keys in the field sections are similar, described in :ref:`Common Field Settings <sec:common_settings>`.
Some specific field keys are shown below:

.. _tab:velocityparams:

.. csv-table:: ``FLUID VELOCITY`` settings in the ``.par`` file
   :widths: 20,20,60
   :header: Key, Value(s), Description/Note(s)/Default Value
  
   ``density`` / ``rho``,``<float>``, "Fluid density"
   ``viscosity`` / ``mu``,``<float>``, "Fluid dynamic viscosity"


.. _tab:scalarparams:

.. csv-table:: ``SCALAR FOO`` settings in the ``.par`` file (specific to scalar ``FOO``)
   :widths: 20,20,60
   :header: Key, Value(s), Description/Note(s)/Default Value
  
   ``mesh``,``fluid`` |br| ``+ solid``, "Specifies the mesh region where scalar ``FOO`` is solved (relevant to :term:`CHT` case) |br| Default = ``fluid``"
   ``transportCoeff``,``<float>``, "Transport property for the scalar ``FOO`` (e.g., :math:`\rho c_p` for ``TEMPERATURE``) in the ``fluid`` ``mesh``"
   ``diffusionCoeff``,``<float>``, "Diffusion coefficient for the scalar ``FOO`` (e.g., :math:`k` for ``TEMPERATURE``) in the ``fluid`` ``mesh``"
   ``transportCoeffSolid``,``<float>``, "Transport property for the scalar ``FOO`` (e.g., :math:`\rho c_p` for ``TEMPERATURE``) in the ``solid`` ``mesh``"
   ``diffusionCoeffSolid``,``<float>``, "Diffusion coefficient for the scalar ``FOO`` (e.g., :math:`k` for ``TEMPERATURE``) in the ``solid`` ``mesh``"

.. _sec:common_settings:

Common Field Settings
^^^^^^^^^^^^^^^^^^^^^

The following table describes settings and corresponding keys for the linear solver.
The keys are common to all solution fields, including velocity, pressure and scalar fields.
These are to be included in the ``.par`` file under appropriate section for ``FLUID VELOCITY``, ``FLUID PRESSURE``, general ``SCALAR`` and specific scalar (``SCALAR FOO``).

.. note::

   Linear solver settings for all scalar fields can be commonly specified under the ``SCALAR`` section.
   Any setting under the specific ``SCALAR FOO`` section will override the common settings under ``SCALAR`` for ``FOO`` field

.. _tab:commonparams:

.. csv-table:: Common settings for all fields in the ``.par`` file
   :widths: 20,20,60
   :header: Key, Value(s), Description/Note(s)/Default Value

   ``solver``,"``none`` |br| ``user`` |br| ``cvode`` |br| ``CG`` |br| ``+ combined`` |br| ``+ block`` |br| ``+ flexible`` |br| ``+ maxiter=<int>`` |br| ``GMRES`` |br| ``+ flexible`` |br| ``+ maxiter=<int>`` |br| ``+ nVector=<int>`` |br|  ``+ iR``","Solve off |br| user-specified |br| CVODE solver (see :ref:`sec:cvodepars`) |br| Conjugate gradient solver. **Default solver for velocity and scalar equation** |br| **Default for scalar equation** |br| **Default velocity solver** |br| . |br| . |br| . |br| Generalized Minimal Residual solver. **Default solver for pressure** |br| **Default for pressure** |br| . |br| Dimension of Krylov space |br| Iterative refinment "  
   ``residualTol``,"``<float>`` |br| ``+ relative=<float>``","absolute linear solver residual tolerance. Default = ``1e-4`` |br| use absolute/relative residual (whatever is reached first)"
   ``absoluteTol``,"``<float>``","absolute solver tolerance (for CVODE only) |br| Default = ``1e-6``"
   ``initialGuess``,"``previous`` |br| ``extrapolation`` |br| ``projection`` |br| ``projectionAconj`` |br| ``+ nVector=<int>``", ". |br| **Default for velocity and scalars** |br| . |br| Defaults for pressure |br| dimension of projection space"
   ``preconditioner``,"``Jacobi`` |br| ``multigrid`` |br| ``+ multiplicative`` |br| ``+ additive`` |br| ``+ SEMFEM`` |br| ``SEMFEM``","**Default for velocity and scalars** |br| Polynomial multigrid + coarse grid projection. **Default for pressure** |br| Default |br| . |br| smoothed SEMFEM |br| ."
   ``coarseGridDiscretization``,"``FEM`` |br| ``+ Galerkin`` |br| ``SEMFEM``","Linear finite element discretization. Default |br| coarse grid matrix by Galerkin projection |br| Linear FEM approx on high-order nodes"
   ``coarseSolver/semfemSolver``,"``smoother`` |br| ``jpcg`` |br| ``+ residualTol=<float>`` |br| ``+ maxiter=<int>`` |br| ``boomerAMG`` |br| ``+ smoother`` |br| ``+ cpu`` |br| ``+ device`` |br| ``+ overlap``", ". |br| Jacobi preconditioned CG |br| . |br| . |br| Hypre's AMG solver |br| . |br| . |br| . |br| overlap coarse grid solve in additive MG cycle"
   ``pMGSchedule``,"``p=<int>, degree=<int>, ...``","custom polynomial order and Chebyshev order for each pMG level"
   ``smootherType``,"``Jacobi`` |br| ``ASM, RAS`` |br| ``+ Chebyshev`` |br| ``+ FourthChebyshev`` |br| ``+ FourthOptChebyshev`` |br| ``+ maxEigenvalueBoundFactor=<float>``",". |br| overlapping additive/restrictive Schwarz |br| 1st Kind Chebyshev acceleration |br| 4th Kind Chebyshev acceleration |br| 4th Opt Chebyshev acceleration |br| ."
   ``checkPointing``, ``true``/``false``, "Turns on/off checkpointing for specific field |br| Default = ``true``"
   ``boundaryTypeMap``,"``<bcType for ID 1>, <bcType for ID 1>, ...``","See :ref:`boundary_conditions` for details"
   ``regularization``,"``hpfrt`` |br| ``+ nModes=<int>`` |br| ``+ scalingCoeff=<float>`` |br| ``gjp`` |br| ``+ scalingCoeff=<float>`` |br| ``avm`` |br| ``+ c0`` |br| ``+ scalingCoeff=<float>`` |br| ``+ noiseThreshold=<float>`` |br| ``+ decayThreshold=<float>`` |br| ``+ activationWidth=<float>``","High-pass filter stabilization |br| number of modes |br| filter strength |br| Gradient Jump Penalty |br| scaling factor in penalty factor fit |br| Artificial Viscosity Method |br| make viscosity C0 |br| . |br| smaller values will be considered to be noise |br| . |br| half-width of activation function"

.. _sec:cvodepars:

CVODE Parameters
"""""""""""""""""""""
.. _tab:cvodeparams:

.. csv-table:: ``CVODE`` settings in the ``.par`` file
   :widths: 20,20,60
   :header: Key, Value(s), Description/Note(s)/Default Value

   ``solver``,"``cbGMRES, GMRES`` |br| ``+ nVector=<int>``", "Linear solver |br| Dimension of Krylov space"
   ``gsType``,"``classical, modified``", ""
   ``relativeTol``,"``<float>``", "relative tolerance |br| Default = ``1e-4``"
   ``epsLin``,``<float>``,"ratio between linear and nonlinear tolerances |br| Default = ``0.5``"
   ``dqSigma``,``<float>``,"step size for Jv difference quotient |br| Default = ``automatic``"
   ``maxSteps``,``<int>``,""
   ``sharedRho``,"``true`` / ``false``", "use same *density* field for all but the first scalar |br| Default = ``false``"
   ``jtvRecycleProperties``,"``true`` / ``false``","recycle property (freeze) evaluation for Jv |br| Default = ``true``"
   ``dealiasing``,"``true`` / ``false``",""

.. _sec:neknekpars:

NekNek Parameters
"""""""""""""""""""""
.. _tab:neknekparams:

.. csv-table:: ``NEKNEK`` settings in the ``.par`` file
   :widths: 20,20,60
   :header: Key, Value(s), Description/Note(s)/Default Value

   ``boundaryEXTOrder``,``<int>``, "Boundary extrapolation order |br| Default = ``1``. >1 may require additional corrector steps"
   ``multirateTimeStepping``,"``true, false`` |br| ``+ correctorSteps=<int>``","Default = ``false`` |br| Outer corrector steps. Default is ``0``. Note: ``boundaryEXTOrder`` > 1 requires ``correctorSteps`` > 0 for stability"
   

.. _udf_functions:

User-Defined Host File (.udf)
-----------------------------

The ``.udf`` file is a :term:`OKL` and C++ mixed language source file, where user code used to formulate the case is placed.
This code is placed in various user-defined functions (*UDFs*) and these can be used to perform virtually any action that can be programmed in C++.
Some of the more common examples are setting initial conditions, querying the solution at regular intervals, and defining custom material properties and source terms.
The available functions that you may define in the ``.udf`` file are as follows.

.. _okl_block:

OKL block
"""""""""

The ``.udf`` typically includes a ``#ifdef __okl__`` block which is where all OKL code is placed that runs on the compute backend specified to :term:`OCCA`.
The most frequent use of this block is to provide the functions for boundary conditions that require additional information, such as a value to impose for a Dirichlet velocity condition, or a flux to impose for a Neumann condition.
Additional user functions may be placed in this block to allow advanced modification of the simulation or post-processing functionality, such as calculating exact values at a specified time point.
Example generic skeleton of typical code structure in :term:`OKL` block is shown below:

.. code-block::
  
  #ifdef __okl__

  @kernel void computeexact(const dlong Ntotal)
  {
    for (dlong n = 0; n < Ntotal; ++n; @tile(p_blockSize, @outer, @inner)) {
      if (n < Ntotal) {
        // some code
      }
    }
  }

  void udfDirichlet(bcData \*bc)
  {
    if(isField("fluid velocity")) {
      bc->uxFluid = 1.0;
      bc->uyFluid = 0.0;
      bc->uzFluid = 0.0;
    }
    else if (isField("fluid pressure")) {
      bc->pFluid = 0.0;
    }
    else if (isField("scalar temperature")) {
      bc->sScalar = 0.0;
    }
  }

  void udfNeumann(bcData \*bc)
  {
    if(isField("fluid velocity")) {
      bc->tr1 = 0.0;
      bc->tr2 = 0.0;
    }
    else if (isField("scalar temperature")) {
      bc->fluxScalar = 0.0;
    }
  }

  #endif

.. tip::

  If the user-defined functions are sufficiently large, it is conventional practice to write them in a ``.oudf`` file which is included within the ``ifdef`` block instead of the functions in the ``.udf`` file, as follows:

  .. code-block:: c++

     #ifdef __okl__

     #include "case.oudf"

     #endif

Details of the ``udfDirichlet`` and ``udfNeumann`` functions used for setting Dirichlet and Neumann boundary conditions, respectively, can be found in :ref:`boundary_conditions`.

.. _udf_setup0:

UDF_Setup0
""""""""""

This user-defined function is passed the nekRS :term:`MPI` communicator ``comm`` and a data structure containing all of the user-specified simulation options, ``options``.
This function is called once at the beginning of the simulation *before* initializing the nekRS internals such as the mesh, solvers, and solution data arrays.
Because virtually no aspects of the nekRS simulation have been initialized at the point when this function is called, this function is primarily used to modify or read the user settings.
Example usage is show below:

.. code-block:: c++
  
   static dfloat P_GAMMA;

   void UDF_Setup0(MPI_Comm comm, setupAide &options)
   {
     platform->par->extract("casedata","p_gamma",P_GAMMA);
   }

   void UDF_LoadKernels(deviceKernelProperties& kernelInfo)
   {
     kernelInfo.define("p_GAMMA") = P_GAMMA;
   }

In the above example ``UDF_Setup0`` routine is used to read ``p_gamma`` key value defined in user section ``CASEDATA`` in the ``.par`` file (see :ref:`user section <sec:user_section>`).
``platform->par->extract`` is a convenient function available in *NekRS* to perform this operation.
The extracted value is assigned to a global variable ``P_GAMMA`` defined at the top of ``.udf`` file and later assigned to a preprocessor macro, made available on the device kernels during JIT compilation.

UDF_LoadKernels
"""""""""""""""

As shown in the example above, ``UDF_LoadKernels`` is primarily used in the ``.udf`` file to append preprocessor macros (global directives) to kernel files.
It takes an argument ``deviceKernelProperties& kernelInfo`` which stores the metadata for kernel compilation.
``kernelInfo.define`` function is used to define the kernel macros and these can later be used in any of the kernel functions.

UDF_Setup
"""""""""

The ``UDF_Setup`` function is called once at the beginning of the simulation *after* initializing the mesh arrays, solution arrays, material property arrays, and boundary field mappings. 
It is typically the function in ``.udf`` the user will interact with the most. 
Various operations are performed within this routine, including, but not limited to:

* Assign initial conditions (see :ref:`initial_conditions`).
* Mesh manipulation (see :ref:`tutorial_rans` tutorial).
* Assign function pointers to user-defined spatially varying material properties (see :ref:`properties`).
* Assign function pointers to user-defined source terms (see :ref:`source_terms`).
* Initialize and setup RANS turbulence models (see :ref:`ktau_model`)
* Initialize and setup Low-Mach compressible model (see :ref:`lowmach_model`)
* Initialize solution recyling routines and arrays (see :ref:`recycling`)
* Allocate ``bc->o_usrwrk`` array for assigning user-defined boundary conditions (see *TODO*)
* Initialize time averaging routines and arrays (see *TODO*)


UDF_ExecuteStep
"""""""""""""""

This user-defined function provides the most flexibility of all the *NekRS* user-defined functions.
It is called once at the start of the simulation just before beginning the time stepping, and then once per time step after running each step.
Two arguments are passed to this routine, including current time (``double``) and timestep (``int``).
Various operations are performed within this routine, including, but not limited to:

* Call time averaging routines (see *TODO*).
* Call solution recycling routines (see :ref:`recycling`).
* Various post-processing operations:

  * Extracting data over a line (see :ref:`extract_line`).
  * Write custom field files (see *TODO*)

.. _usr_functions:

Legacy Nek5000 User File (.usr)
--------------------------------

*NekRS* provides an optional framework for legacy interface with the *Nek5000* code, allowing access to fortran 77 based *Nek5000* user routines to perform custom operations.
The user has the option to include ``<case>.usr`` in the case directory to include the usual *Nek5000* user routines. 
For users unfamiliar with *Nek5000* code, more information can be found in `Nek5000 documentation <https://nek5000.github.io/NekDoc/>`_.
Note that not all *Nek5000* routines are called by *NekRS*. 
More commonly, the user may require call to the ``userchk()`` routine in *Nek5000* for post-processing operations. 
If required, it must be explicitly called from ``.udf`` file as shown below:

.. code-block:: c++

   void UDF_ExecuteStep(double time, int tstep)
   {
     if(nrs->checkpointStep) {
        nrs->copyToNek(time, tstep);
        nek::userchk();
     }
   }

For most applications, the ``userchk`` routine will be called from ``UDF_ExecuteStep`` function, likely for post-processing operations.
``nrs->copyToNek`` copies all solution fields from :term:`OCCA` arrays to *Nek5000* (fortran) arrays. 
This call is necessary before calling ``nek::userchk`` in order for the user to perform any post-processing on field arrays in *Nek5000*.

.. warning::

   The ``nrs->copyToNek`` call performs expensive operation of copying the data from :term:`OCCA` arrays to Nek5000.
   This must be done sparingly, only at certain time steps in the simulation.
   Remember to call this routine within suitable ``if`` condition block.
   As shown in the example above, ``nrs->copyToNek`` and ``nek::userchk`` are called only at ``checkPointStep``.

Details on *Nek5000* ``.usr`` file can be found `here <https://nek5000.github.io/NekDoc/problem_setup/usr_file.html#user-routines-file-usr>`_ and specific information on ``userchk`` fortran routine `here <https://nek5000.github.io/NekDoc/problem_setup/usr_file.html#userchk>`_.

Other *Nek5000* user routines that are internally called by *NekRS* during initialization are ``usrdat0``, ``usrdat``, ``usrdat2`` and ``usrdat3``.
Details on these initialization routines can be found `here <https://nek5000.github.io/NekDoc/problem_setup/usr_file.html#initialization-routines>`_.
These routines can be optionally used for specifying boundary conditions, mesh manipulation, parameter specification or other initialization operations. 

Legacy Data Interface
"""""""""""""""""""""

*NekRS* provides an in-built mechanism to pass variables or array pointers to share data between *Nek5000* and *NekRS* through ``nekrs_registerPtr`` fortran routine.
Consider the following code snippet in ``.usr`` file:

.. code-block:: fortran

   subroutine userchk()
   include 'SIZE'
   include 'TOTAL'

   common /exact/ uexact(lx1,ly1,lz1,lelt*3),
  &               texact(lx1,ly1,lz1,lelt) 

   real uexact, texact

   call computeexact(uexact, texact)

   call nekrs_registerPtr('uexact', uexact)
   call nekrs_registerPtr('texact', texact)

   return
   end

   subroutine computeexact(uexact, texact)
   include 'SIZE'
   include 'TOTAL'

   ! Code to compute exact solution

   return
   end

   subroutine usrdat0

   real gamma 
   save gamma

   gamma = 1.4

   call nekrs_registerPtr('gamma', gamma)

   return
   end

In the above code, two routines are defined in the fortran common block ``exact`` to store exact solution for velocity and temperature.
The exact solutions are computed in ``userchk`` and the array pointers for the solutions are registered using ``nekrs_registerPtr`` subroutine.
It takes two arguments - a string identifier for the pointer and the pointer to the array to be registered. 
(Note that pointer to the first memory location in the array is registered).
It is critical that these arrays are declared in fortran common block or that they are saved for them to be visible globally.
Another example is shown with a variable, ``gamma``, in ``usrdat0`` routine, which is also registered in a similar manner.
Note that the variable ``gamma`` is made static (or saved in memory) using the ``save`` command. 

These pointers can now be accessed in ``.udf`` file and used to transfer data between *NekRS* and *Nek5000*.
Example usage in ``.udf`` is shown below:

.. code-block:: c++

   static dfloat gamma;

   void UDF_Setup0(MPI_Comm comm, setupAide &options)
   {
      gamma = *nek::ptr<double>("gamma"); 
   }

   void UDF_LoadKernels(deviceKernelProperties& kernelInfo)
   {
      kernelInfo.define("p_GAMMA") = gamma;
   }

   void UDF_ExecuteStep(double time, int tstep)
   {
      if(nrs->lastStep) {
        auto mesh = nrs->meshV;

        nek::userchk();

        std::vector<double> uexact(nek::ptr<double>("uexact"), nek::ptr<double>("uexact") + nrs->fluid->fieldOffsetSum);
        std::vector<double> texact(nek::ptr<double>("texact"), nek::ptr<double>("texact") + mesh->Nlocal);

        auto o_uexact = platform->device.malloc<dfloat>(nrs->fluid->fieldOffsetSum);
        auto o_texact = platform->device.malloc<dfloat>(nrs->fluid->fieldOffset);

        o_uexact.copyFrom(uexact);
        o_texact.copyFrom(texact);

        //compute error here...
      }
   }

As shown, ``nek::ptr`` stores the registered pointers, which is recognised using the string identifier specified in ``.usr`` file.  
In ``UDF_Setup0`` the value referenced by the pointer corresponding to ``gamma`` is assigned to the C++ static variable of the same name.
It is later used to define a kernel macro in ``UDF_LoadKernels`` as shown.
Similarly, the pointers to fortran arrays identified by ``uexact`` and ``texact`` are used to copy data onto ``std::vector`` containers of the same name.
The arrays are then copied from ``std::vector`` containers to :term:`OCCA` arrays ``o_uexact`` and ``o_texact``.
The user can then perform any required operations on these arrays, such as compute solution error norms.

Mesh File (.re2)
----------------

*TODO*

The nekRS mesh file is provided in a binary format with a nekRS-specific
``.re2`` extension. This format can be produced by either:

* Converting a mesh made with commercial meshing software to ``.re2`` format, or
* Directly creating an ``.re2``-format mesh with nekRS-specific scripts

There are three main limitations for the nekRS mesh:

* nekRS is restricted to 3-D hexahedral meshes.
* The numeric IDs for the mesh boundaries must be ordered contiguously beginning from 1.
* The ``.re2`` format only supports HEX8 and HEX 20 (eight- and twenty-node) hexahedral elements.

Lower-dimensional problems can be accommodated on these 3-D meshes by applying zero gradient
boundary conditions to all solution variables in directions perpendicular to the
simulation plane or line, respectively. All source terms and material properties in the
governing equations must therefore also be fixed in the off-interest directions.

For cases with conjugate heat transfer, nekRS uses an archaic process
for differentiating between fluid and solid regions. Rather than block-restricting variables to
particular regions of the same mesh, nekRS retains two independent mesh representations
for the same problem. One of these meshes represents the flow domain, while the other
represents the heat transfer domain. The ``nrs_t`` struct, which encapsulates all of
the nekRS simulation data related to the flow solution, represents the flow mesh as
``nrs_t.mesh``. Similarly,
the ``cds_t`` struct, which encapsulates all of the nekRS simulation data related to the
convection-diffusion passive scalar solution, has one mesh for each passive scalar. That is,
``cds_t.mesh[0]`` is the mesh for the first passive scalar, ``cds_t.mesh[1]`` is the mesh
for the second passive scalar, and so on.
Note that only the temperature passive scalar uses the conjugate heat transfer mesh,
even though the ``cds_t`` struct encapsulates information related to all other
passive scalars (such as chemical concentration, or turbulent kinetic energy). All
non-temperature scalars are only solved on the flow mesh.

.. warning::

  When writing user-defined functions that rely on mesh information (such as boundary
  IDs and spatial coordinates), you must take care to use the correct mesh representation
  for your problem. For instance, to apply initial conditions to a flow variable, you
  would need to loop over the number of quadrature points known on the ``nrs_t`` meshes,
  rather than the ``cds_t`` meshes for the passive scalars (unless the meshes are the same,
  such as if you have heat transfer in a fluid-only domain).
  Also note that the ``cds_t * cds`` object will not exist if your problem
  does not have any passive scalars.

nekRS requires that the flow mesh be a subset of the heat transfer mesh. In other words,
the flow mesh always has less than (or equal to, for cases without conjugate heat transfer)
the number of elements in the heat transfer mesh. Creating a mesh for conjugate heat
transfer problems requires additional pre-processing steps that are described in the
:ref:`Creating a Mesh for Conjugate Heat Tranfser <cht_mesh>` section. The remainder
of this section describes how to generate a mesh in ``.re2`` format, assuming
any pre-processing steps have been done for the special cases of conjugate heat transfer.

.. _trigger_file:

Trigger Files (.upd)
--------------------

**TODO** Full description

Allows modifications to the simulation during execution.
Can be edited and then notify of changes through sending a signal MPI rank 0.
