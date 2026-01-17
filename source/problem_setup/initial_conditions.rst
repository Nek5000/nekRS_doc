.. _initial_conditions:

Initial conditions
==================

There are several locations users can setup the initial conditions. In execution
order, *nekRS* calls ``useric`` in the ``.usr`` file (if a ``.usr`` exists),
reads a checkpoint file (if specified in the ``.par``), calls ``userchk`` (also
only if ``.usr`` exists), and then calls ``UDF_Setup`` in the ``.udf`` file.
After these callbacks, *nekRS* calls ``mesh->update`` to recompute the geometric
factors from the current mesh coordinates. If a checkpoint is used, mesh
coordinates and solution fields from the checkpoint overwrite anything set in
``useric``.

You may still control which fields are read via par-keys and modify coordinates
in ``UDF_Setup``. See :ref:`ic_restart` below for details.
As a final step, the solutions are projected onto a :math:`C^0` space to enforce
continuity across element boundaries.

.. note::

   Everything in ``.usr`` is :ref:`usr_file` in order to reuse *Nek5000*
   settings. If no ``.usr`` file is provided, both ``useric`` and ``userchk``
   are skipped.

.. _ic_udf_setup:

Specification in ``.udf`` File
-------------------------------

``UDF_Setup()`` in :ref:`udf_file` file is the primary location for setting
initial conditons. The minimal example below sets initial conditions for all
three velocity components, temperature and two passive scalars.

.. code-block::

    void UDF_Setup()
    {
      auto mesh = nrs->meshV;

      if (platform->options.getArgs("RESTART FILE NAME").empty()) {
        std::vector<dfloat> U(nrs->fluid->fieldOffsetSum, 0.0);
        std::vector<dfloat> temp(mesh->Nlocal, 0.0);
        std::vector<dfloat> ps1(mesh->Nlocal, 0.0);
        std::vector<dfloat> ps2(mesh->Nlocal, 0.0);

        auto [x, y, z] = mesh->xyzHost();

        for (int n = 0; n < mesh->Nlocal; n++) {
          U[n + 0 * nrs->fieldOffset] = sin(x[n]) * cos(y[n]) * cos(z[n]);
          U[n + 1 * nrs->fieldOffset] = -cos(x[n]) * sin(y[n]) * cos(z[n]);
          U[n + 2 * nrs->fieldOffset] = 0.0;
          temp[n] = x[n];
          ps1[n] = y[n] + 0.5;
          ps2[n] = 1.0;
        }
        nrs->fluid->o_U.copyFrom(U.data(), U.size());
        nrs->scalar->o_solution("temperature").copyFrom(temp.data(), temp.size());
        nrs->scalar->o_solution("ps1").copyFrom(ps1.data(), ps1.size());
        nrs->scalar->o_solution("ps2").copyFrom(ps2.data(), ps2.size());
      }
    }

Because ``UDF_Setup`` is called after checkpoint file(s) are read, the block
above is guarded by ``platform->options.getArgs("RESTART FILE NAME").empty()``
so it only executes when no restart file is used.

``U, temp, ps1, ps2`` are temporary host arrays stored in ``std::vector``
containers. The call ``xyzHost()`` host view of the ``x,y,z`` coordinate, which
are then used to defube spatially varying initial conditions.
Finally, these host arrays are copied into the corresponding device arrays,
``nrs->fluid->o_U`` and ``nrs->scalar->o_solution`` via ``copyFrom`` so that
the initial conditions are available to the simulation.

The scalar identifiers ``"temperature"``, ``"ps1"`` and ``"ps2"`` must match
the names declared in the :ref:`General Section <sec:generalpars>` of the
``.par`` file (see :ref:`Parameter file section for details <par_file>` for
details).

.. tip::

   You can also overwrite solution fields in ``UDF_ExecuteStep``, for example
   to impose a time-dependent analytic solution and initialize all lagged
   stages during the first few timesteps. See the ``ethier`` example for
   details.

.. _ic_restart:

Restarting from Field File(s)
-------------------------------

The default restart mechanism is quite capable and supports:

1. Simulation and checkpoint files using different polynomial orders.
2. Restarts with a different number of MPI ranks.
3. Built-in :ref:`h-refinement <meshing_hrefine>`, allowing restarts from a
   previously unrefined mesh.
4. Selecting fields from multiple checkpoint files.
5. Grid-to-grid interpolation, so checkpoint files may come from an entirely
   different mesh provided with the coordinates are stored in the file.

Items 1-3 are handled automatically. Options 4 and 5 can be specified as
restart options, controlled by the par-key ``startFrom`` (see
:ref:`tab:generalparams`).

.. csv-table:: Examples of restart options
   :widths: 60,40
   :header: "``startFrom``", "Description"
   :class: tall

   "``""ethier0.f00000""``", "Read all fields from the file, including time."
   "``""ethier0.f00000""+U+time=2.0``", "Read only the velocity and set time to 2.0."
   "``""ethier0.f00000""+X, ""ethier0.f00001""+U``", "Read the mesh from one file and velocity from another."
   "``""ethier0.f00000+int""``", "Interpolate all fields from a checkpoint written on a different mesh."
   "``""ethier.bp""``", "Read all fields from the last step in the bp file."
   "``""ethier.bp""+U+step=0``", "Read the velocity from the first step in bp file."
   "``""ethier.bp""+int``", "Interpolate all fields from a bp written on different mesh."

.. note::

   ADIOS2 bp files may store multiple time steps in the same folder. By default,
   *nekRS* reads the last available step, but you can select a specific one via
   the ``step=<int>`` option in ``startFrom``.

.. tip::

   *nekRS* supports two checkpoint formats, and both are binary. However, you
   can still inspect the state of a file.

   - Nek-format checkpoint file (``#.f#####``):
     See the `Nek5000 documentation <https://nek5000.github.io/NekDoc/problem_setup/case_files.html#restart-output-files-f>`_ for details.
     You can examine the first 132 ASCII characters of its header. For example:

     .. code-block:: bash

        head -c 132 ethier0.f00000

   - ADIOS2 bp5 format:
     You can inspect the metadata with the ``bpls`` executable. See
     :ref:`postproc_adios` for usage.

Users can also manually load checkpoint files in the ``.udf`` by calling
``nrs->restartFromFiles`` inside functions such as ``UDF_Setup`` or
``UDF_ExecuteStep`` in the following format.

.. code-block:: cpp

   std::vector<std::string> fileList = {"ethier0.f00000+X", "ethier0.f00001+U+P+S"}
   nrs->restartFromFiles(fileList);

.. _ic_usr_file:

Legacy ``useric`` in ``.usr`` File
----------------------------------

See the `Nek5000 documentation <https://nek5000.github.io/NekDoc/problem_setup/usr_file.html#useric>`__
for example usage.

