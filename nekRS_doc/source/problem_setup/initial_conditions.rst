.. _initial_conditions:

Initial conditions
==================

The initial conditions are specified in ``UDF_Setup()`` routine.
The following minimal example code snippet sets initial conditions for all three components of velocity, temperature and two passive scalars.

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
    
The above code block is not executed if the simulation is restarted from a prior field file, in which case the initial conditions are read from that file.
It is essential to perform initialization only if ``platform->options.getArgs("RESTART FILE NAME").empty()`` is true, otherwise the restart field information will be over-written.
``U, temp, ps1, ps2`` are temporary array containers declared on host memory to assign initial conditions. 
The call ``xyzHost()`` returns the ``x,y,z`` coordinate arrays which may be used to assign spatially varying initial conditions as shown above.
Finally, the host arrays must be copied into device arrays, ``nrs->fluid->o_U`` and ``nrs->scalar->o_solution`` using ``copyFrom`` calls for initial conditions to be made available for subsequent simulation.
Note that the scalar identifiers ``"temperature"``, ``"ps1"`` and ``"ps2"`` must be identical to the string declared in ``.par`` file (see :ref:`Parameter file section for details <parameter_file>`)

