.. _boundary_conditions:

Boundary conditions
===================

Boundary conditions for each mesh boundary should normally be set in the :ref:`parameter_file`, using the ``boundaryTypeMap`` parameter.
This is used within the ``FLUID VELOCITY``, ``SCALAR TEMPERATURE`` or ``SCALAR XXXXX`` sections to set the boundary conditions of the respective solvers.

Available Types
---------------

All available boundary conditions are summarised in the command line help function for nekRS (``nekrs --help par``) and below.

.. literalinclude:: ../../parHelp.txt
   :language: none
   :lines: 1-3, 154-182

.. tip::

   To setup some cases (e.g., cases that use third party meshes wherein boundaries are identified by a unique integer), you may need to use the ``boundaryIDMap`` parameter of the ``MESH`` section to apply the ``boundaryTypeMap`` options to the correct boundary IDs of the mesh (by default NekRS assumes that the boundaryIDs start at 1).
   Assuming an inlet-outlet pipe flow, below is an example where the three boundary conditions are applied to the boundary IDs 389 (``udfDirichlet`` - inlet), 231 (``zeroNeumann`` - outlet), 4 (``zeroDirichlet`` - walls).

   .. code-block::

      [MESH]
      boundaryIDMap = 389, 231, 4

      [FLUID VELOCITY]
      boundaryTypeMap = udfDirichlet, zeroNeumann, zeroDirichlet


.. note::

  For conjugate heat transfer cases, it is necessary to also specify the velocity boundary map, ``boundaryIDMapFluid``, separately to identify the interface conformal walls between the fluid and solid mesh.
  Following the above example, assume that the pipe is enclosed in a solid annular insulated jacket where the boundary ID 4 corresponds to the conformal wall between the solid and fluid mesh, while remaining insulated walls of the solid mesh are 10, 20 and 30.

    .. code-block::

       [MESH]
       boundaryIDMapFluid = udfDirichlet, zeroNeumann, zeroDirichlet
       boundaryIDMap = zeroNeumann, zeroNeumann, zeroNeumann, udfDirichlet, zeroNeumann

       [FLUID VELOCITY]
       boundaryTypeMap = 389, 231 ,4

       [SCALAR TEMPERATURE]
       boundaryTypeMap = 10, 20, 30, 389, 231

User Defined Value/Gradient
"""""""""""""""""""""""""""

If a boundary condition requires specifying a value rather than just a type (specified by the ``boundaryTypeMap`` entries ``udfDirichlet`` or ``udfNeumann`` in ``.par`` file), the value must be defined in the ``udfDirichlet`` or ``udfNeumann`` function calls within the ``.oudf`` file (or :ref:`okl block <okl_block>` section of ``.udf`` file), for Dirichlet or Neumann boundary types, respectively.
``bcData`` struct is passed as an argument to these functions, which has the following parameters available:

+-----------------------------------------------------+------------------------------+---------------------------------------+
| Name                                                | Type                         | Description                           |
+=====================================================+==============================+=======================================+
| ``fieldOffset``                                     | ``int``                      | array offset                          |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``idxVol``                                          | ``int``                      | vol index of boundary GLL             |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``id``                                              | ``int``                      | boundaryID                            |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``time``                                            | ``double``                   | current simulation time               |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``x``, ``y``, ``z``                                 | ``dfloat``                   | location coordinates of boundary GLL  |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``nx``, ``ny``, ``nz``                              | ``dfloat``                   | local normal vector components        |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``t1x``, ``t1y``, ``t1z``                           | ``dfloat``                   | local tangent vector components       |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``t2x``, ``t2y``, ``t2z``                           | ``dfloat``                   | local bi-tangent vector components    |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``tr1``                                             | ``dfloat``                   | local traction along tangent          |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``tr2``                                             | ``dfloat``                   | local traction along bi-tangent       |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``uxFluid``, ``uyFluid``, ``uzFluid``               | ``dfloat``                   | local velocity components             |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``pFluid``                                          | ``dfloat``                   | local pressure                        |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``uxFluidInt``, ``uyFluidInt``, ``uzFluidInt``      | ``dfloat``                   | interpolated velocity components      |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``idScalar``                                        | ``int``                      | scalar index                          |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``sScalar``                                         | ``dfloat``                   | scalar value (Dirichlet)              |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``fluxScalar``                                      | ``dfloat``                   | scalar flux (Neumann)                 |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``sScalarInt``                                      | ``dfloat``                   | interpolated scalar value             |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``sInfScalar``                                      | ``dfloat``                   | ?                                     |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``h``                                               | ``dfloat``                   | ?                                     |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``uxGeom``, ``uyGeom``, ``uzGeom``                  | ``dfloat``                   | mesh velocity components              |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``transCoeff``                                      | ``dfloat``                   | :math:`\rho` or :math:`\rho c_p`      |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``diffCoeff``                                       | ``dfloat``                   | diffusion coefficient                 |
+-----------------------------------------------------+------------------------------+---------------------------------------+
| ``usrwrk``                                          | ``@globalPtr const dfloat*`` | pointer to user array                 |
+-----------------------------------------------------+------------------------------+---------------------------------------+

To specify the user defined boundary condition to the corresponding field, the field is identified in the udf functions using the ``isField()`` parameter.
The ``isField()`` macro is defined in the ``bcData`` struct and it takes field string identifier as the parameter. 
The string identifiers are:

+-----------------------------------+-------------------------------------------+
| String Identifier                 | Field                                     |
+===================================+===========================================+
| ``"fluid velocity"``              | fluid                                     |
+-----------------------------------+-------------------------------------------+
| ``"scalar temperature"``          | temperature                               |
+-----------------------------------+-------------------------------------------+
| ``"scalar XXX"``                  | field corresponding to scalar ``"XXX"``   |
+-----------------------------------+-------------------------------------------+

.. note::

    The scalar string identifiers are declared in :ref:`Parameter file <parameter_file>`.


.. tip::

    Some example cases are listed below to show the usage of different user defined boundary condition types:

    +------------------+-------------------------------------+------------------------------------------------------------------------------------+
    | BC type          | Field                               | Example case                                                                       |
    +------------------+-------------------------------------+------------------------------------------------------------------------------------+
    | ``udfDirichlet`` | ``"fluid velocity"``                | `ethier <https://github.com/Nek5000/nekRS/tree/next/examples/ethier>`_             |
    |                  +-------------------------------------+------------------------------------------------------------------------------------+
    |                  | ``"fluid velocity"`` (recycling)    | `turbPipe <https://github.com/Nek5000/nekRS/tree/next/examples/turbPipe>`_         |
    |                  +-------------------------------------+------------------------------------------------------------------------------------+
    |                  | ``"fluid velocity"`` (interpolated) | `eddyNekNek <https://github.com/Nek5000/nekRS/tree/next/examples/eddyNekNek>`_     |
    |                  +-------------------------------------+------------------------------------------------------------------------------------+
    |                  | ``"fluid pressure"``                | `hemi <https://github.com/Nek5000/nekRS/tree/next/examples/hemi>`_                 |
    |                  +-------------------------------------+------------------------------------------------------------------------------------+
    |                  | ``"scalar XXX"``                    | `conj_ht <https://github.com/Nek5000/nekRS/tree/next/examples/conj_ht>`_           |
    |                  +-------------------------------------+------------------------------------------------------------------------------------+
    |                  | ``"scalar XXX"`` (interpolated)     | `eddyNekNek <https://github.com/Nek5000/nekRS/tree/next/examples/eddyNekNek>`_     |
    |                  +-------------------------------------+------------------------------------------------------------------------------------+
    |                  | ``"geom"``                          | `mv_cyl <https://github.com/Nek5000/nekRS/tree/next/examples/mv_cyl>`_             |
    +------------------+-------------------------------------+------------------------------------------------------------------------------------+
    | ``udfNeumann``   | ``"fluid velocity"``                | `gabls1 <https://github.com/Nek5000/nekRS/tree/next/examples/gabls1>`_             |
    |                  +-------------------------------------+------------------------------------------------------------------------------------+
    |                  | ``"scalar XXX"``                    | `ethier <https://github.com/Nek5000/nekRS/tree/next/examples/ethier>`_             |
    +------------------+-------------------------------------+------------------------------------------------------------------------------------+    

Internal / Periodic
"""""""""""""""""""
.. _periodic_boundary:

``None`` is used when a internal boundary condition is required or a periodic boundary condition has been set as part of the mesh and it does not need to be considered as part of the standard processing of boundary conditions.
Perodicity is linked to the mesh connectivity and is handled by the meshing tool.


Turbulent Inflow (for LES/DNS)
------------------------------

Velocity Recycling
""""""""""""""""""

NekRS offers an in-built technique to impose fully developed turbulent inflow conditions.  
It comprises recycling the velocity from an offset surface in the domain interior to the inlet boundary resulting in a quasi-periodic domain. 
The required routines are available through the ``planarCopy.hpp`` header file. 
Sample ``.udf`` code to use velocity recycling is as follows,

.. code-block::

  #include "planarCopy.hpp"

  deviceMemory<int> o_bID;
  planarCopy* recyc = nullptr;
  dfloat area, uBulk;

  #ifdef __okl__
   void udfDirichlet(bcData * bc)
   {
    if (isField("fluid velocity")) {
      bc->uxFluid = bc->usrwrk[bc->idxVol + 0 * bc->fieldOffset];
      bc->uyFluid = bc->usrwrk[bc->idxVol + 1 * bc->fieldOffset];
      bc->uzFluid = bc->usrwrk[bc->idxVol + 2 * bc->fieldOffset];
    }
   }
  #endif

  void UDF_Setup()
  {
    auto mesh = nrs->meshV;
    const auto zmax = platform->linAlg->max(mesh->Nlocal, mesh->o_z, platform->comm.mpiComm());
    const auto zmin = platform->linAlg->min(mesh->Nlocal, mesh->o_z, platform->comm.mpiComm());
    const auto zLen = abs(zmax - zmin);
                    
    platform->app->bc->o_usrwrk.resize(mesh->dim * nrs->fluid->fieldOffset);

    uBulk = 1.0;                            // mean bulk velocity
    const dfloat Dh = 1.0;                  // hydraulic diameter
    const dfloat xOffset = 0.0;             // x-offset of recycling surface 
    const dfloat yOffset = 0.0;             // y-offset of recycling surface
    const dfloat zOffset = 10.0 * Dh;       // z-offset of recycling surface

    std::vector<int> bID;
    bID.push_back(1);
    o_bID.resize(bID.size());
    o_bID.copyFrom(bID);

    deviceMemory<dfloat> o_tmp(mesh->Nlocal);
    platform->linAlg->fill(o_tmp.size(), 1.0, o_tmp);
    area = mesh->surfaceAreaMultiplyIntegrate(o_bID, o_tmp);

    recyc = new planarCopy(mesh, 
                           nrs->fluid->o_solution(),
                           mesh->dim,
                           nrs->fluid->fieldOffset,
                           xOffset,
                           yOffset,
                           zOffset,
                           bID[0],
                           platform->app->bc->o_usrwrk);
  }

  void UDF_ExecuteStep(double time, int tstep)
  {
    recyc->execute();

    const auto flux = mesh->surfaceAreaNormalMultiplyVectorIntegrate(nrs->fluid->fieldOffset, o_bID, platform->app->bc->o_usrwrk);
    platform->linAlg->scale(nrs->fluid->fieldOffsetSum, -uBulk * area / flux, platform->app->bc->o_usrwrk);
  }

The above example assumes that the inlet surface, with boundary ID (``bID``) equal to 1, is perpendicular to the *x-y* plane.
The recycling surface is located at an offset distance equal to :math:`20 Dh` from the inlet surface, where :math:`Dh` is the hydraulic diameter of the inlet cross-section. 
``uBulk`` is the desired mean velocity at the inlet boundary. 

.. note::

  The recommended recycling length is :math:`8-10 Dh` [Lund1998]_ to allow the resolution of largest eddies. The user must ensure that the boundary layer is not acted upon by significant pressure gradients or geometrical changes between the inlet and recycle plane.

The interpolated velocities from the recycling surface are stored in the ``o_usrwrk`` array which is accessible in the ``udfDirichlet`` kernel.
Therefore, ``o_usrwrk`` vector must be resized to accomodate the velocity data to ``mesh->dim * nrs->fluid->fieldOffset`` length.
The ``planarCopy()`` call setups up the necessary interpolation procedures to map the velocity from domain interior to the inlet surface.
Note the parameters passed as arguments to the setup routine.
The ``execute()`` routine is called from the ``UDF_ExecuteStep`` routine at every time step.
It, internally, populates the ``o_usrwrk`` array with the scaled instantaneous velocity at the recycle plane based on specified ``ubulk``.
Finally, the captured velocity is specified as inlet Dirichlet boundary condition in ``udfDirichlet`` in the :ref:`okl block <okl_block>` section, as shown above.
Note the proper indexing and offset specified using ``bc->idxVol`` and ``bc->fieldOffset`` for the three velocity components.
