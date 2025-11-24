.. _other_resources:

Other Resources
===============

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
   function, and :math:`F` a differential operator. The strong form appears in
   integrals of the type

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

Table :numref:`tab:opsem` lists the ``opSEM`` functions and their corresponding
discrete operators. The following notation is used:

- :math:`D_x, D_y, D_z`: first-order derivative operators in the
  :math:`x`, :math:`y`, and :math:`z` directions.
- :math:`B`: multiplication by the (lumped) mass matrix.
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
   auto o_gradU = opSEM::strongGradVec(mesh, nrs->fieldOffset, nrs->fluid->o_U);

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

