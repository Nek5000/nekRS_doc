.. _running:

Running
=======

This page gives information on how to run *NekRS* after it has been installed (see :ref:`Installation instructions <installing>` ) and appropriate input files have been generated 
(see :ref:`case`).

Native Running
--------------

To run *NekRS* natively from the directory containing the :ref:`case`, use the following general command

.. code-block::

    mpirun -np <number of MPI tasks> nekrs --setup <par|sess file> [options]

Below are the command line options/argument that can be used to further modifiy how *NekRS* is run.

.. _tab:runoptions:

.. csv-table:: Command line options to run *NekRS*
   :widths: 20,40,40,10
   :header: Command line argument, Options, Description,Required

   ``--help``,"None |br| ``par`` |br| ``env``", "Prints summary of available command line arguments |br| Prints summary of :ref:`Parameter file <par_file>` sections and keys |br| Prints list of *NekRS* enviorenment variables","No"
   ``--setup``,``par`` file |br| ``sess`` file, "Name of ``.par`` file of the case |br| Name of ``.sess`` file (NekNek)","Yes"
   ``--build-only``,"None |br| ``#procs``","Runs JIT compilation for ``np`` procs |br| Runs JIT compilation for specified number of processors","No"
   ``--cimode``,"``<id>``","Runs *NekRS* in CI mode corresponding to specified ``<id>``","No"
   ``--debug``,"None", Run *NekRS* in debug mode,"No"
   ``--attach``,"None",?,"No"
   ``--backend``,"``CPU`` |br| ``CUDA`` |br| ``HIP`` |br| ``DPCPP`` |br| ``OPENCL``","Manually set the backend device for running","No"
   ``--output``,"``<local-path-to-logfile>``",Name of log file to print *NekRS* output,"No"
   ``--device-id``,"``<id>`` |br| ``LOCAL-RANK``",Manually set OCCA device ID for GPU |br| For CPU,"No"

.. _nekrs_scripts:

MPI launch scripts
^^^^^^^^^^^^^^^^^^^

A number of scripts ship with *NekRS* itself and are located in the ``$NEKRS_HOME/bin``. A brief summary of these scripts and their usage is as follows.

* ``nrsmpi <casename> <processes>``: run *NekRS* in parallel with ``<processes>`` parallel processes for the case files that are prefixed with ``casename``.
* ``nrsbmpi <casename> <processes>``: same as ``nrsmpi``, except that nekRS runs in the background and the output is printed to log file

Running on HPC Systems
----------------------

Specific scripts and instructions for running on specific HPC systems can be found `in nekRS_HPCsupport Github repository <https://github.com/Nek5000/nekRS_HPCsupport>`_
