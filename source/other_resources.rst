.. _other_resources:

Other Resources
===============

.. _occa_memory:

OCCA Memory
-----------

*TODO*

``deviceMemory``
``occa::memory``
memory Pool

.. _mesh_t:

Mesh Class (``mesh_t``)
-----------------------

Each solver in *NekRS* owns a corresponding ``mesh_t`` object.

- For the fluid (velocity / pressure) solver, the mesh is stored in
  ``nrs->fluid->mesh`` (often also referred to as the fluid mesh).
- For scalar fields (e.g., temperature), meshes are accessed through the
  scalar, for example,
  ``nrs->scalar->mesh("temperature")`` or by index via
  ``nrs->scalar->mesh(idx)``, where the ``idx`` can be accessed via
  ``idx = nrs->scalar->nameToIndex.find("temperature")``

In most cases, all scalar meshes point to the same fluid-domain mesh. However,
in CHT setups, selected scalars (such as temperature)
can live on a mesh that covers both fluid and solid regions.
(See :ref:`conjugate heat transfer turorial <conjugate_heat_transfer>`)

*NekRS* performs a domain decomposition to distribute the mesh across MPI ranks
with approximately equal degrees of freedom, in order to balance the
computational load. As a result, the members of ``mesh_t`` correspond to the
**local** portion of the mesh assigned to a given MPI rank. For example,
``Nelements`` is the number of elements owned locally by that rank.

In the summary table below, we refer generically to the ``mesh`` object
(e.g., ``mesh->Nelements``) without distinguishing whether it is the fluid
mesh, a scalar mesh, or a CHT mesh. The meaning of each member is the same,
but applied to the local partition held by that particular solver instance.

If there are both host and device, the device one will have prefix ``o_``

.. table:: Important ``mesh_t`` members
   :name:  mesh_data

   =========================== ================================== ================== =================================================
   Variable Name               Size                               Host or Device?           Meaning
   =========================== ================================== ================== =================================================
   ``dim = 3``                 1                                  Host               spatial dimension of mesh
   ``Nfaces = 6``              1                                  Host               number of faces per element
   ``Nverts = 8``              1                                  Host               number of vertices per element
   ``NfaceVertices = 4``       1                                  Host               number of vertices per face
   ``N``                       1                                  Host               polynomial order for each dimension
   ``Nq = N + 1``              1                                  Host               number of quadrature points in each direction
   ``Np = Nq * Nq * Nq``       1                                  Host               number of quadrature points per element
   ``Nfp = Nq * Nq``           1                                  Host               number of quadrature points per face
   ``Nfields``                 1                                  Host               number of fields passed to the PDE solver
   ``Nelements``               1                                  Host               number of local elements owned by current process
   ``Nlocal = Nelements * Np`` 1                                  Host               number of local quadrature points owned by current process
   ``NboundaryFaces``          1                                  Host               *total* number of faces on a boundary (rank sum)

   ``solid``                   1                                  Host               whether contains solid domain (1 = CHT)
   ``elementInfo``             ``Nelements``                      Both               phase of element (0 = fluid, 1 = solid)

   ``EToV``                    ``Nelements * Nverts``             Host               element to vertex connectivity
   ``EToE``                    ``Nelements * Nfaces``             Host               element to element connectivity
   ``EToF``                    ``Nelements * Nfaces``             Host               element to local face connectivity
   ``EToP``                    ``Nelements * Nfaces``             Host               element to partition/process connectivity
   ``EToB``                    ``Nelements * Nfaces``             Both               mapping of elements to type of boundary condition
   ``vmapM``                   ``Nelements * Nfaces * Nfp``       Both               quadrature point index for faces on boundaries

   ``ogs``                     ``ogs_t``                          ogs_t              OCCA gather/scatter operation handle
   ``oogs``                    ``oogs_t``                         oogs_t             (overlapped) OCCA gather/scatter operation handle

   ``o_x``                     ``Nelements * Np``                 device             :math:`x`-coordinates of physical quadrature points
   ``o_y``                     ``Nelements * Np``                 device             :math:`y`-coordinates of physical quadrature points
   ``o_z``                     ``Nelements * Np``                 device             :math:`z`-coordinates of physical quadrature points

   ``gllz``                    ``Nq``                             Both               1D GLL points
   ``gllw``                    ``Nq``                             Both               1D GLL quadratule weights
   ``D``                       ``Nq * Nq``                        Both               1D Differentiation matrix

   ``o_LMM (or o_Jw)``         ``Nlocal``                         Device             lumped mass matrix
   ``o_invLMM (or o_invAJw)``  ``Nlocal``                         Device             inverse of assembled lumped mass matrix
   ``o_vgeo``                  ``Nlocal * 12``                    Device             volume geometric factors
   ``o_sgeo``                  ``Nelements * Nfaces * Nfp * 19``  Device             surface geometric factors
   ``o_ggeo``                  ``Nlocal * 7``                     Device             second-order volumetric geometric factors
   =========================== ================================== ================== =================================================

.. note::

   For the mapping from reference coordinates :math:`(r,s,t)` to physical
   coordinates :math:`(x,y,z)`, the geometric factors are stored in three main
   device arrays:

   - ``o_vgeo``: volume geometric factors

     - the 9 components :math:`\partial (r,s,t) / \partial (x,y,z)`,
     - the Jacobian :math:`J`,
     - the weighted Jacobian :math:`J w`,
     - and its inverse :math:`1 / (J w)`.

   - ``o_sgeo``: surface geometric factors (on face nodes).
     Only 18 of the 19 stored entries are currently used:

     - the 3 components of the outward unit normal,
     - 6 components for the tangential and bitangential vectors,
     - the surface Jacobian and weighted surface Jacobian,
     - 3 entries for :math:`n \cdot \nabla r`, :math:`n \cdot \nabla s`,
       :math:`n \cdot \nabla t`,
     - 3 entries for averaged normal/tangential data.

   - ``o_ggeo``: geometric factors for the 3D Laplacian

     - the 6 unique components of the symmetric metric tensor :math:`G_{ij}`,
     - and the weighted Jacobian :math:`J w`.

.. _linalg:

Linear Algebra Functions (linAlg)
---------------------------------

doxygen?

*TODO*

``src/platform/linAlg/linAlg.tpp``
``src/platform/linAlg/linAlg.hpp``
``src/platform/linAlg/linAlg.cpp``



.. _opsem:

SEM Operators (opSEM)
---------------------

The ``opSEM`` namespace (see ``src/core/opSEM.hpp`` for arguments) collects
commonly used spectral-element derivative-based operators. By default, the
operators correspond to the weak form, while functions with a ``strong`` suffix
implement the corresponding weighted strong form. Strong functions come
with an optional argument ``avg`` (default: ``true``) that returns the
unweighted version.

.. note::

   In the finite element setting, let :math:`u` be a field, :math:`v` a test
   function, and :math:`F` a differential operator. The (weighted) strong form
   appears in integrals of the type

   .. math::

      \int_\Omega v\, F[u] \,\mathrm{d}\vec{x}.

   If :math:`F` can be written in divergence form
   :math:`F[u] = \nabla \cdot \widehat{F}[u]`, with the associated flux
   :math:`\widehat{F}[u]`, integration by parts gives

   .. math::

      \int_\Omega v\, F[u] \,\mathrm{d}\vec{x}
      = - \int_\Omega (\nabla v) \cdot \widehat{F}[u] \,\mathrm{d}\vec{x}
        + \int_{\partial\Omega} v\, \widehat{F}[u]\cdot\vec{n}\,\mathrm{d}s.

   Neglecting the boundary term, the weak form is the volume term with a minus
   sign. In ``opSEM``, the weak operators implement

   .. math::

      \int_\Omega (\nabla v) \cdot \widehat{F}[u] \,\mathrm{d}\vec{x},

   i.e., the same expression without the minus sign.

   After discretization, testing with all basis functions turns this bilinear
   form into a linear system for the coefficients (nodal values) of :math:`u`.
   The test functions :math:`v` no longer appear explicitly in the algebraic
   system.

:numref:`tab:opsem` lists the ``opSEM`` functions and their corresponding
discrete operators. The following notations are used:

- :math:`D_x, D_y, D_z`: first-order derivative operators in the
  :math:`x`, :math:`y`, and :math:`z` directions.
- :math:`B`: multiplication by the diagonal (lumped) mass matrix.
- :math:`G`: gather–scatter averaging,
  :math:`G = M^{-1} Q Q^{T}`, where :math:`Q` is the finite-element gather
  matrix, :math:`Q^{T}` is the scatter, and :math:`M` is a diagonal matrix of
  nodal multiplicities.
- :math:`\Lambda`: multiplication by a scalar field :math:`\lambda`
  (a diagonal operator in the discrete setting).

.. _tab:opsem:

.. csv-table:: ``opSEM`` functions applied to a discretized scalar :math:`u`
               or a vector :math:`\vec{u} = (u_x, u_y, u_z)`
   :widths: 20,10,10,30,30
   :class: tall
   :header: Function,In,Out,Math,Avg mode (``avg=true``)

   ``grad``, "Scalar", "Vector (3)", ":math:`[B D_x u,\; B D_y u,\; B D_z u]`", "NA"
   ``strongGrad``, "Scalar", "Vector (3)", "Same as ``grad``", ":math:`[G D_x u,\; G D_y u,\; G D_z u]`"
   ``strongGradVec``, "Vector (3)", "Vector (9)", "``strongGrad`` of :math:`u_x`, then of :math:`u_y` and of :math:`u_z`", "Applying ``strongGrad(avg=true)`` to 3 components"
   ``divergence``, "Vector (3)", "Scalar", ":math:`B (D_x^{T} u_x + D_y^{T} u_y + D_z^{T} u_z)`", "NA"
   ``strongDivergence``, "Vector (3)", "Scalar", ":math:`B (D_x u_x + D_y u_y + D_z u_z)`", ":math:`G (D_x u_x + D_y u_y + D_z u_z)`"
   ``laplacian``, "Scalar", "Scalar", ":math:`(D_x^{T} \Lambda B D_x + D_y^{T} \Lambda B D_y + D_z^{T} \Lambda B D_z)\, u`", "NA"
   ``strongLaplacian``, "Scalar", "Scalar", ":math:`B \big(D_x \Lambda G D_x + D_y \Lambda G D_y + D_z \Lambda G D_z\big) u`", ":math:`G \big(D_x \Lambda G D_x + D_y \Lambda G D_y + D_z \Lambda G D_z\big) u`"
   ``strongCurl``, "Vector (3)", "Vector (3)", ":math:`[B (D_y u_z - D_z u_y),\; B (D_z u_x - D_x u_z),\; B (D_x u_y - D_y u_x)]`", ":math:`[G (D_y u_z - D_z u_y),\; G (D_z u_x - D_x u_z),\; G (D_x u_y - D_y u_x)]`"

The derivative operators commonly used in postprocessing are listed in
:numref:`tab:opsem_lookup`.

.. _tab:opsem_lookup:

.. csv-table:: Lookup table for ``opSEM`` functions
   :widths: 40,60
   :class: tall
   :header: Operators, ``opSEM`` function

   ":math:`\nabla u`",                        "``strongGrad(avg=true)``"
   ":math:`\nabla \vec u`",                   "``strongGradVec(avg=true)``"
   ":math:`\nabla \cdot \vec u`",             "``strongDivergence(avg=true)``"
   ":math:`\nabla \cdot (\lambda \nabla u)`", "``strongLaplacian(avg=true)``"
   ":math:`\nabla\times\vec u`",              "``strongCurl(avg=true)``"

For example, the following call computes the gradient of the velocity field:

.. code-block:: cpp

   auto mesh = nrs->meshV;
   auto o_gradU = opSEM::strongGradVec(mesh, nrs->fieldOffset, nrs->fluid->o_U); // default: avg = true

The output ``o_gradU`` is a device buffer from the memory pool, with type
``deviceMemory<dfloat> o_out(mesh->dim * mesh->dim * nrs->fieldOffset);``
and its entries ordered as

.. math::

   \partial_x u_x,\; \partial_y u_x,\; \partial_z u_x,\;
   \partial_x u_y,\; \partial_y u_y,\; \partial_z u_y,\;
   \partial_x u_z,\; \partial_y u_z,\; \partial_z u_z.




.. _nek5000_tools:

Building the Nek5000 Tools
--------------------------

*NekRS* does not package the legacy toolkit available with :term:`Nek5000`, e.g., tools for creating or adapting meshes.
It relies instead on the toolbox available with :term:`Nek5000`, which must be installed separately.
To build these scripts, clone the Nek5000 repository, and then navigate to the ``tools`` directory.
The desired tools can be built here using `./maketools`.

For example, if you want to make the ``genbox`` tool,

.. code-block:: bash

  cd ~
  git clone https://github.com/Nek5000/Nek5000.git
  cd Nek5000/tools
  ./maketools genbox

This will create binary executables in the ``Nek5000/bin`` directory.
You may want to add this folder to your ``PATH`` environment variable for quick access to these tools,

.. code-block:: bash

    export PATH+=:$HOME/Nek5000/bin

Additional information about *Nek5000* tools can be found `here <https://nek5000.github.io/NekDoc/tools.html>`_.

.. note::

  The `genbox <https://nek5000.github.io/NekDoc/tools/genbox.html>`_ tool creates simple 2D and 3D meshes that can be used for both *Nek5000* and *NekRS* simulations.
  For *NekRS* applications, this is most likely the only tool that be of use to most users.
  Other *Nek5000* tools will be rarely used.

.. note::

   Some tools require additional dependencies such as ``cmake``, ``wget`` (and
   internet access), or X11 libraries. Check the messages from the
   ``maketools`` script and the corresponding ``build.log`` files in each tool
   directory for any errors.

