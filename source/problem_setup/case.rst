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
* **Session file**, with ``.sess`` extension. This file is for NekNek only.
  (See :ref:`overlapping_overset_grids` tutorial).

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

  # This is a comment
  [FOO]
    key = value # this is also a comment
    baz = bat

  [BAR]
    alpha = beta
    gamma = delta + keyword=value + ...   # composed value with inline keywords
    theta = value1, value2, value3        # comma-separated list

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

  - Section name and key/value pairs are case insensitive
  - Values are words with all spaces removed
  - Values enclosed within quotes preserve case and whitespace
  - Values prefixed with 'env::' are interpreted as references to environment variables
  - Values separated by commas forms a list

.. _sec:user_section:

User Sections
""""""""""""""""""

Custom sections may be added via the ``userSections`` key to pass additional
keys into the code.

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
   :class: tall
   :header: Key, Value(s), Description/Note(s)/Default Value

   ``polynomialOrder``,``<int>``, "``polynomialOrder`` > 10 is currently not supported"
   ``dealiasing``,``true`` / ``false``, "Enables/disables over-integration of convective term |br| Default = ``true``"
   ``cubaturePolynomialOrder``,``<int>``, "Polynomial order of ``dealiasing`` |br| Default = 3/2*(``polynomialOrder`` + 1) - 1"
   ``verbose``,``true`` / ``false``, "``true`` instructs *NekRS* to print detailed diagnostics to *logfile* |br| Default = ``false``"
   ``redirectOutputTo``,``<string>``,"String entry for the name of the *logfile* to direct *NekRS* output"
   ``startFrom``,"``""<string>"", ...`` |br| ``+ time=<float>`` |br| ``+ x`` |br| ``+ u`` |br| ``+ s or s00 s01 s02 ...`` |br| ``+ int``", "Restart from specified ``<string>`` file |br| reset ``time`` to specified value |br| read mesh coordinates |br| read velocity |br| read all scalar or specified scalars |br| interpolate solution (useful if mesh coordinates are different)"
   ``timeStepper``,``tombo1`` / ``tombo2`` / ``tombo3``,"Order of time discretization for BDFk/EXTk scheme |br| Default = ``tombo2``"
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
   ``udf``,"``""<string>""``","Optional name of user-defined host function file |br| Default is ``<case>.udf``"
   ``oudf``,"``""<string>""``","Optional name of user-defined OCCA kernel function file |br| As a default *NekRS* expects these are defined in :ref:`OKL block <okl_block>` in ``.udf`` file"
   ``usr``,"``""<string>""``","Optional name of user-defined legacy *Nek5000* (fortran) function file |br| Default is ``<case>.usr``"
   ``regularization``,"","Specifies regularization options for all fields |br| See :ref:`common field settings<sec:common_settings>` for details"

.. _sec:occa:

OCCA Parameters
""""""""""""""""
.. _tab:occaparams:

.. csv-table:: ``OCCA`` keys in the ``.par`` file
   :widths: 20,20,60
   :class: tall
   :header: Key, Value(s), Description/Note(s)/Default Value

   ``backend``, |br| ``SERIAL`` / |br| ``CUDA`` / |br| ``HIP`` /|br| ``DPCPP``,"Specifies the *device* for JIT compilation. Default is defined in ``$NEKRS_HOME/nekrs.conf`` |br| CPU |br| NVIDIA GPU (CUDA) |br| AMD GPU (HIP) |br| Intel GPU (oneAPI)"
   ``deviceNumber``,``<int>`` |br| ``LOCAL-RANK``,"Default is ``LOCAL-RANK``"
   ``platformNumber``,``<int>``, "Only used by ``DPCPP`` |br| Default is ``0``"

.. _sec:problemtype:

Problem Type Parameters
""""""""""""""""""""""""""
.. _tab:problemparams:

.. csv-table:: ``PROBLEMTYPE`` keys in the ``.par`` file
   :widths: 20,20,60
   :class: tall
   :header: Key, Value(s), Description/Note(s)/Default Value

   ``equation``,``stokes`` |br| ``navierStokes`` |br| ``+ variableViscosity``, "Stokes solver |br| Navier-Stokes solver |br| uses stress formulation (required for spatially varying viscosity)"

.. _sec:meshpars:

Mesh Parameters
""""""""""""""""
.. _tab:meshparams:

.. csv-table:: ``MESH`` keys in the ``.par`` file
   :widths: 20,20,60
   :class: tall
   :header: Key, Value(s), Description/Note(s)/Default Value

   ``partitioner``,``rbc`` / ``rsb`` / ``rbc+rsb``,"Specifies mesh partitioner |br| Default = ``rbc+rsb`` "
   ``boundaryIDMap``,"``<int>, <int>, ...``", "Map mesh boundary ids to 1,2,3,... |br| See :ref:`boundary conditions<boundary_conditions>` for details"
   ``boundaryIDMapFluid``,"``<int>, <int>, ...``", "Required for conjugate heat transfer cases |br| See :ref:`boundary conditions<boundary_conditions>` for details"
   ``connectivityTol``,"``<float>``","Specifies mesh tolerance for partitioner |br| Default = ``0.2``"
   ``file``,"``""<string>""``","Optional name of mesh (``.re2``) file |br| Default is ``<case>.re2``"


.. _sec:field_settings:

Field Settings
"""""""""""""""""""""

The sections for specific fields, including velocity (``FLUID VELOCITY``), pressure (``FLUID PRESSURE``) and scalars (``SCALAR`` or ``SCALAR FOO``) contain keys to describe linear solver setting for the corresponding field.
Most of the keys in the field sections are similar, described in :ref:`Common Field Settings <sec:common_settings>`.
Some specific field keys are shown below:

.. _tab:velocityparams:

.. csv-table:: ``FLUID VELOCITY`` settings in the ``.par`` file
   :widths: 20,20,60
   :class: tall
   :header: Key, Value(s), Description/Note(s)/Default Value

   ``density`` / ``rho``,``<float>``, "Fluid density"
   ``viscosity`` / ``mu``,``<float>``, "Fluid dynamic viscosity"


.. _tab:scalarparams:

.. csv-table:: ``SCALAR FOO`` settings in the ``.par`` file (specific to scalar ``FOO``)
   :widths: 20,20,60
   :class: tall
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
   :widths: 15,35,50
   :class: tall
   :header: Key, Value(s), Description/Note(s)/Default Value

   ``solver``,"``none`` |br| ``user`` |br| ``cvode`` |br| ``CG`` |br| ``+ combined`` |br| ``+ block`` |br| ``+ flexible`` |br| ``+ maxiter=<int>`` |br| ``GMRES`` |br| ``+ flexible`` |br| ``+ maxiter=<int>`` |br| ``+ nVector=<int>`` |br|  ``+ iR``","Solve off |br| user-specified |br| CVODE solver (see :ref:`sec:cvodepars`) |br| Conjugate gradient solver. **Default solver for velocity and scalar equation** |br| **Default for scalar equation** |br| **Default velocity solver** |br| . |br| . |br| . |br| Generalized Minimal Residual solver. **Default solver for pressure** |br| **Default for pressure** |br| . |br| Dimension of Krylov space |br| Iterative refinment"
   ``residualTol``,"``<float>`` |br| ``+ relative=<float>``","absolute linear solver residual tolerance. Default = ``1e-4`` |br| use absolute/relative residual (whatever is reached first)"
   ``absoluteTol``,"``<float>``","absolute solver tolerance (for CVODE only) |br| Default = ``1e-6``"
   ``initialGuess``,"``previous`` |br| ``extrapolation`` |br| ``projection`` |br| ``projectionAconj`` |br| ``+ nVector=<int>``", ". |br| **Default for velocity and scalars** |br| . |br| Defaults for pressure |br| dimension of projection space"
   ``preconditioner``,"``Jacobi`` |br| ``multigrid`` |br| ``+ multiplicative`` |br| ``+ additive`` |br| ``+ SEMFEM`` |br| ``SEMFEM``","**Default for velocity and scalars** |br| Polynomial multigrid + coarse grid projection. **Default for pressure** |br| Default |br| . |br| smoothed SEMFEM |br| ."
   ``coarseGridDiscretization``,"``FEM`` |br| ``+ Galerkin`` |br| ``SEMFEM``","Linear finite element discretization. Default |br| coarse grid matrix by Galerkin projection |br| Linear FEM approx on high-order nodes"
   ``coarseSolver/semfemSolver``,"``smoother`` |br| ``jpcg`` |br| ``+ residualTol=<float>`` |br| ``+ maxiter=<int>`` |br| ``boomerAMG`` |br| ``+ smoother`` |br| ``+ cpu`` |br| ``+ device`` |br| ``+ overlap``", ". |br| Jacobi preconditioned CG |br| . |br| . |br| Hypre's AMG solver |br| . |br| . |br| . |br| overlap coarse grid solve in additive MG cycle"
   ``pMGSchedule``,"``p=<int> + degree=<int>, ...``","custom polynomial order and Chebyshev order for each pMG level"
   ``smootherType``,"``Jacobi`` |br| ``ASM, RAS`` |br| ``+ Chebyshev`` |br| ``+ FourthChebyshev`` |br| ``+ FourthOptChebyshev`` |br| ``+ maxEigenvalueBoundFactor=<float>``",". |br| overlapping additive/restrictive Schwarz |br| 1st Kind Chebyshev acceleration |br| 4th Kind Chebyshev acceleration |br| 4th Opt Chebyshev acceleration |br| ."
   ``checkPointing``, ``true``/``false``, "Turns on/off checkpointing for specific field |br| Default = ``true``"
   ``boundaryTypeMap``,"``<bcType for ID 1>, <bcType for ID 2>, ...``","See :ref:`boundary_conditions` for details"
   ``regularization``,"``hpfrt`` |br| ``+ nModes=<int>`` |br| ``+ scalingCoeff=<float>`` |br| ``gjp`` |br| ``+ scalingCoeff=<float>`` |br| ``avm`` |br| ``+ c0`` |br| ``+ scalingCoeff=<float>`` |br| ``+ noiseThreshold=<float>`` |br| ``+ decayThreshold=<float>`` |br| ``+ activationWidth=<float>``","High-pass filter stabilization |br| number of modes |br| filter strength |br| Gradient Jump Penalty |br| scaling factor in penalty factor fit |br| Artificial Viscosity Method |br| make viscosity C0 |br| . |br| smaller values will be considered to be noise |br| . |br| half-width of activation function"

.. _sec:cvodepars:

CVODE Parameters
"""""""""""""""""""""
.. _tab:cvodeparams:

.. csv-table:: ``CVODE`` settings in the ``.par`` file
   :widths: 20,20,60
   :class: tall
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
   :class: tall
   :header: Key, Value(s), Description/Note(s)/Default Value

   ``boundaryEXTOrder``,``<int>``, "Boundary extrapolation order |br| Default = ``1``. >1 may require additional corrector steps"
   ``multirateTimeStepping``,"``true, false`` |br| ``+ correctorSteps=<int>``","Default = ``false`` |br| Outer corrector steps. Default is ``0``. Note: ``boundaryEXTOrder`` > 1 requires ``correctorSteps`` > 0 for stability"


.. _udf_functions:

User-Defined Host File (.udf)
-----------------------------

A *UDF* (user-defined function) is an entry point NekRS calls during setup and
time stepping, letting a case customize behavior without touching core code.
Typical uses include setting initial and boundary conditions, defining custom
materials, models, and source terms, and performing sampling and logging for
post-processing. The available entry-point functions and their typical uses are
as follows. All UDF functions are **optional**. If a function is not provided,
NekRS uses a default no-op version.

.. _okl_block:

OKL block
"""""""""

The ``.udf`` usually contains a preprocessor guard ``#ifdef __okl__`` that
encloses all :term:`OKL` kernels and device functions compiled for the selected
:term:`OCCA` backend. For an overview of the OKL language, see the `OKL language guide
<https://libocca.org/#/okl/introduction>`_.
NekRS provides default device functions for boundary conditions.
Details of the ``udfDirichlet`` and ``udfNeumann`` functions used for Dirichlet
and Neumann boundary conditions, respectively, can be found in :ref:`boundary_conditions`.

.. code-block:: cpp

  #ifdef __okl__

  void udfDirichlet(bcData *bc)
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

  void udfNeumann(bcData *bc)
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


Functions marked with the decorate ``@kernel`` are compiled as device kernels
and can be launched from UDF host code directly as a regular function.

.. code-block:: cpp

   // -------- Device side (inside the OKL block) --------
   #ifdef __okl__

   @kernel void my_exact(const dlong Ntotal,
                         @ restrict const dfloat *X,
                         @ restrict dfloat *U)
   {
     for (dlong n = 0; n < Ntotal; ++n; @tile(p_blockSize, @outer, @inner)) {
       if (n < Ntotal) {
         U[n] = sin(X[n]);
       }
     }
   }

   #endif

   // -------- Host side (UDF code) --------
   {
     auto mesh = nrs->meshV;
     deviceMemory<dfloat> o_tmp(mesh->Nlocal);
     my_exact(o_tmp.size(), mesh->o_x, o_tmp); // o_tmp = sin(x)
   }

.. tip::

   If the user-defined functions are sufficiently large, it is conventional practice to write them in a ``.oudf`` file which is included within the ``ifdef`` block instead of the functions in the ``.udf`` file, as follows:

  .. code-block:: c++

     #ifdef __okl__

     #include "case.oudf"

     #endif

.. tip::

   Many common operations are available in ``platform->linAlg`` and ``opSEM``
   (e.g., fills, axpby, dot products, norms, gradient, divergence). Prefer these
   utilities over writing custom kernels when possible. See :ref:`linalg` and
   :ref:`opsem` for details.

.. _udf_setup0:

UDF_Setup0
""""""""""

``UDF_Setup0`` is called **once** at startup, before *NekRS* initializes the
mesh, solvers, or solution arrays. It receives the *NekRS* :term:`MPI`
communicator (``comm``) and the user options object (``options``). Because the
simulation state is not yet constructed, this function is mainly for reading or
adjusting user settings.

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

In this example, ``UDF_Setup0`` reads the ``p_gamma`` key from the ``CASEDATA``
user section in the ``.par`` file (see :ref:`sec:user_section`) using the
convenience helper ``platform->par->extract``. The value is stored in a global
``P_GAMMA`` and then exported as the device macro ``p_GAMMA`` in
``UDF_LoadKernels`` for JIT-compiled kernels.

.. warning::
   Do **not** modify parameters used to compile legacy *Nek5000* in
   ``UDF_Setup0`` such as ``polynomialOrder``.


UDF_LoadKernels
"""""""""""""""

As shown in the example above, ``UDF_LoadKernels`` is primarily used to attach
preprocessor macros (global defines) for device kernels.
It takes ``deviceKernelProperties& kernelInfo`` which holds compilation metadata.
Using ``kernelInfo.define("NAME") = value`` to expose constants to :term:`OKL`.

.. tip::

   For constants passed to device code (via ``UDF_LoadKernels``), it’s recommended
   to use a ``p_`` prefix (e.g., ``p_gamma``, ``p_Pr``, ``p_Re``) for clarity and
   consistency.

UDF_Setup
"""""""""

``UDF_Setup`` is called once at the beginning of the simulation *after* *nekRS*
initializing the mesh, solution fields, material properties, and boundary mappings.
It is typically the ``.udf`` function users will use most to overwrite defaults
and customize settings.
Various operations are performed within this routine, including, but not limited to:

* Assign initial conditions (see :ref:`initial_conditions`).
* Mesh manipulation (see :ref:`tutorial_rans` tutorial).
* Register function pointers for user-defined spatially varying material properties (see :ref:`properties`).
* Register function pointers for user-defined source terms (see :ref:`source_terms`).
* Initialize and set up RANS turbulence models (see :ref:`ktau_model`)
* Initialize and set up the Low-Mach compressible model (see :ref:`lowmach_model`)
* Initialize solution recycling routines and arrays (see :ref:`recycling`)
* Allocate ``bc->o_usrwrk`` for user-defined boundary data (see *TODO*)
* Initialize time-averaging routines and buffers (see *TODO*)


UDF_ExecuteStep
"""""""""""""""

``UDF_ExecuteStep`` offers the most flexibility. It is called once just before
time marching begins, and then once **per time step**. The routine receives the
current time (``double t``) and the step index (``int tstep``).
Typical operations include:

* Run time-averaging updates (see *TODO*)
* Invoke solution recycling logic (see :ref:`recycling`).
* Post-processing tasks, such as:

  * Extracting data over a line (see :ref:`extract_line`).
  * Writing custom field files (*TODO*).


.. _usr_functions:

Legacy Nek5000 User File (.usr)
--------------------------------

NekRS includes an optional legacy interface to *Nek5000* so you can reuse
Fortran-77 user routines. If ``<case>.usr`` is present in the case directory,
the file is compiled and linked, though not all *Nek5000* user subroutines are
invoked by NekRS.

Many one-time setup routines are still honored so existing *Nek5000* settings can
be carried over. Users can continually use subroutines
``usrdat0``, ``usrdat``, ``usrdat2``, ``usrdat3``, ``usrsetvert`` and ``useric``
to modify mesh geometry, tag boundary surfaces, adjust connectivity and set initial
conditions.
Compute intensive subroutines such as ``userbc``, ``userf``, ``userq`` and ``uservp``
are **not** used by *NekRS*. Users need to convert them into UDF or OKL implementation.

The ``userchk`` is **not** called automatically. If needed, call it manually
from UDF (e.g., ``UDF_ExecuteStep``) for post-processing or other setup tasks,
for example:

.. code-block:: c++

   void UDF_ExecuteStep(double time, int tstep)
   {
     if(nrs->checkpointStep) { // occasionally call when a checkpoint step
        nrs->copyToNek(time, tstep); // push device field to Nek5000 arrays
        nek::userchk();              // call userchk in .usr
        nrs->copyFromNek(time);      // only if userchk modify solutions
     }
   }

.. note::

   - ``nrs->copyToNek`` copies all solution fields from :term:`OCCA` device arrays
     to *Nek5000* (Fortran) host arrays. Use it only when ``nek::userchk``
     needs current solutions.
   - This transfer is expensive. Wrap the call under a suitable condition to
     avoid excessive overhead.
   - Call ``nrs->copyFromNek`` only if the legacy routine modifies solution
     data that must be synchronized back to the device.

For background on Nek5000, see the `Nek5000 documentation <https://nek5000.github.io/NekDoc/>`_.
Details on these initialization routines are documented `here <https://nek5000.github.io/NekDoc/problem_setup/usr_file.html#initialization-routines>`_.
To migrate a case, follow the *TODO* tutorial to converting a *Nek5000* case to
*NekRS*.
In this section, we focus one the interface for referencing variables between
``.udf`` and ``.usr``.

Legacy Data Interface
"""""""""""""""""""""

NekRS does not “send” values or arrays to *Nek5000*. Instead, it shares them by
**registering pointers**. On the Fortran side, any variable or array you expose
must have a stable address, so put it in a ``COMMON`` block or give it the
``SAVE`` attribute.
Use the Fortran routine ``nekrs_registerPtr`` to register variables/arrays so
NekRS can reference them by a string key.

Here is an example ``.usr`` file:

.. code-block:: fortran

   c-----------------------------------------------------------------------
         subroutine usrdat0
         real gamma
         save gamma

         gamma = 1.4

         call nekrs_registerPtr('gamma', gamma)

         return
         end
   c-----------------------------------------------------------------------
         subroutine userchk()
         include 'SIZE'
         include 'TOTAL'

         common /exact/ uexact(lx1,ly1,lz1,lelt*3),
        &               texact(lx1,ly1,lz1,lelt)

         real uexact, texact

         call nekrs_registerPtr('uexact', uexact(1,1,1,1,1))
         call nekrs_registerPtr('vexact', uexact(1,1,1,1,2))
         call nekrs_registerPtr('wexact', uexact(1,1,1,1,3))
         call nekrs_registerPtr('texact', texact)

         ! update uexact/texact here

         return
         end

In this example, ``COMMON /exact/`` block holds two arrays store exact solutions
(e.g., velocity and temperature). We register the **first-element addresses** of
each velocity components and of temperature with ``nekrs_registerPtr`` inside
``userchk``. The scalar ``gamma`` is registered in ``usrdat0`` and kept alive
with ``SAVE``. Any data you register must be in a ``COMMON`` block or ``SAVE``
so its address remains valid.

These pointers can then be looked up in the ``.udf`` and used to exchange data
between *NekRS* and *Nek5000*. Example usage in ``.udf``:

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
     if (nrs->lastStep) {
       auto mesh = nrs->meshV;
       auto Nlocal = mesh->Nlocal;

       nek::userchk(); // compute exact solutions in .usr

       // Host pointers exported from usr
       auto *u_ex = nek::ptr<double>("uexact");
       auto *v_ex = nek::ptr<double>("vexact");
       auto *w_ex = nek::ptr<double>("wexact");

       // Range constructor (copy values into std::vector)
       std::vector<double> t_ex(nek::ptr<double>("texact"),
                                nek::ptr<double>("texact") + Nlocal);

       // Host buffer
       std::vector<dfloat> uexact(Nlocal), uexact(Nlocal), wexact(Nlocal), texact(Nlocal);

       // Convert from double to dfloat
       for (int i = 0; i < Nlocal; i++) {
         u_exact[i] = static_cast<dfloat>(u_ex[i]);
         v_exact[i] = static_cast<dfloat>(v_ex[i]);
         w_exact[i] = static_cast<dfloat>(w_ex[i]);
         t_exact[i] = static_cast<dfloat>(t_ex[i]);
       }

       // Device buffer
       auto o_uexact = platform->device.malloc<dfloat>(nrs->fluid->fieldOffsetSum);
       auto o_texact = platform->device.malloc<dfloat>(Nlocal);

       // copyFrom(src pointer, count, dest offset, src offset)
       o_uexact.copyFrom(uexact.data(), Nlocal, 0 * nrs->fluid->fieldOffset, 0);
       o_uexact.copyFrom(vexact.data(), Nlocal, 1 * nrs->fluid->fieldOffset, 0);
       o_uexact.copyFrom(wexact.data(), Nlocal, 2 * nrs->fluid->fieldOffset, 0);
       o_texact.copyFrom(texact.data(), Nlocal);

       //compute error here...
     }
   }

As shown, ``nek::ptr`` returns the registered pointers by the string keys
defined in the ``.usr`` file. In ``UDF_Setup0``, the value referenced by
``gamma`` is read into a C++ static and later exported as a device macro in
``UDF_LoadKernels``.

For the arrays (``uexact``, ``vexact``, ``wexact``, ``texact``), we show two ways:
use a raw host pointer from ``nek::ptr<double>(...)`` for velocity components
and use a range constructor to make a copy into ``std::vector`` for temperature.
Then, all are converted from ``double`` to ``dfloat`` before copying to
:term:`OCCA` arrays.

.. note::

   *NekRS* compiles Fortran ``real`` as ``real*8`` (i.e., C++ ``double``).
   In *NekRS*, ``dfloat`` may be ``double`` (default) or ``float`` (FP32 mode).
   Converting on the host as shown works for both settings. If you always run in
   FP64, you can skip the cast. Alternatively, you can convert on device via
   ``platform->copyDoubleToDfloatKernel``.

.. note::

   **Nek5000/NekRS sizes and offsets**

   - ``mesh->Nlocal``: number of local points on this rank
     (e.g., ``lx1*ly1*lz1*nelv``; fluid+solid mesh uses``nelt``).

   - ``fieldOffset``: padded device stride (alignment + max local size).
     It is **not** equal to ``lx1*ly1*lz1*lelv``.

   - ``fieldOffsetSum``: sum of component strides. |br|
     For velocity: ``fieldOffsetSum = 3 * fieldOffset``. |br|
     For scalar, ``fieldOffsetSum = nrs->Nscalar * fieldOffset``.

   Use ``fieldOffset``/``fieldOffsetSum`` for declaring multi-component device
   arrays and indexing. Use ``Nlocal`` for loop lengths and buffer sizes.


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


.. _session_file:

NekNek Session File (.sess)
---------------------------

NekNek allows coupling of multiple overlapping subdomains, relaxing the need for
conformal meshes. Each *nekRS* instance uses its own `.par` file.
In the ``.sess`` file, each line lists one instance as:
``<path-to-par-file>:<num-mpi-ranks>;``

Paths may be absolute or relative, and the sum of ranks over all instances must
equal the total MPI ranks requested. Here is the ``eddyNekNek.sess`` example:

.. code-block:: bash

   inside/inside:1;
   outside/outside:1;

.. tip::

   Using a subfolder per case is optional but recommended. On HPC systems, try
   to size each session in units of a node. For example, if a node has 4 GPUs,
   it’s often best to make each session’s MPI ranks a multiple of 4.


.. _trigger_file:

Trigger Files (.upd)
--------------------

**TODO** Full description

Allows modifications to the simulation during execution.
Can be edited and then notify of changes through sending a signal MPI rank 0.
