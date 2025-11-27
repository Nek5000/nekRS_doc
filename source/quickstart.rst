.. _quickstart:

Quickstart
==========

.. _qstart_before:

-------------------
Before You Begin...
-------------------

This guide assumes basic familiarity with Linux (or another Unix-like OS) and common ``bash`` commands.
The steps below target building *nekRS* on a local machine (desktop or laptop). A **GPU is not required.**

.. note::

   Specific instructions for installing on select :term:`HPC` systems can be
   found `here <https://github.com/Nek5000/nekRS_HPCsupport>`_. We recommend
   becoming familiar with the standard installation first before using these
   HPC-specific instructions.

Before you begin, note the following **requirements** for *NekRS*:

.. mdinclude:: _includes/README.md
   :start-line: 32
   :end-line: 36

Most of these should either be available by default in your OS of choice, or can be installed using a common package manager.

.. tabs::

    .. tab:: Debian/Ubuntu

        Debian based systems (such as Ubuntu) use the ``apt`` package manager.
        GNU Compilers for C++ and fortran alongside compatible OpenMPI and CMake installs can be acquired with:

        .. code-block:: bash

            sudo apt update
            sudo apt install build-essential libopenmpi-dev cmake

    .. tab:: Mac

        The `Homebrew <https://brew.sh/>`_ package manager is commonly used on Mac to provide similar functionality to Linux package managers.
        GNU Compilers for C++ and fortran alongside compatible OpenMPI and CMake installs can be acquired with:

        .. code-block:: bash

            brew install gcc open-mpi cmake

**TODO**: Instructions for installing GPU compilers

.. _installing:

----------------
Installing NekRS
----------------

Acquiring the Code
^^^^^^^^^^^^^^^^^^

*NekRS* can be acquired by either downloading a release:

.. code-block:: bash

    cd ~
    wget https://github.com/Nek5000/nekRS/archive/refs/tags/v25.0-rc1.tar.gz
    tar -xzvf v25.0-rc1.tar.gz
    mv nekRS-25.0-rc1 nekRS
    cd nekRS

or by cloning the GitHub repository:

.. code-block:: bash

   cd ~
   git clone https://github.com/Nek5000/nekRS.git
   cd nekRS
   git checkout next

.. note::

   The latest development version of *NekRS* is maintained on the ``next``
   branch of the GitHub repository. Make sure you have checked out ``next`` (as
   shown above) and verify that your local HEAD matches the latest commit listed at
   `NekRS next branch commits <https://github.com/Nek5000/nekRS/commits/next/>`__.
   You can inspect your local HEAD with ``git log``.

Building NekRS
^^^^^^^^^^^^^^

Use the helper script ``./build.sh`` from the *nekRS* directory.

.. code-block:: bash

   CC=mpicc CXX=mpic++ FC=mpif77 ./build.sh -DCMAKE_INSTALL_PREFIX=$HOME/.local/nekrs

- Adjust C, C++ and Fortran compilers via the environment variables
  ``CC``, ``CXX`` and ``FC``.
- Append any CMake options after ``build.sh``.
  (e.g., ``-DCMAKE_INSTALL_PREFIX`` to specify the installation location.)

.. note::

    The environment variables and CMake options are optional.
    ``CC=mpicc``, ``CXX=mpic++`` and ``FC=mpif77`` are usually the default system ``MPI`` wrapper installed using the ``apt`` or ``brew`` package manager as shown :ref:`above <qstart_before>`.
    If MPI comes from a custom or vendor stack, set these variable accordingly.

.. tip::

   Like standard CMake workflows, *nekrs* configures and compiles in ``./build/`` and installs to ``CMAKE_INSTALL_PREFIX``.
   If you updates *NekRS* or your compiler stacks, make sure to remove the previous ``build/`` and install directory.

The build script firsts configure CMake and prints a summary like:

.. code-block:: none

  ----------------- Summary -----------------
  Installation directory: /home/usr/.local/nekrs
  plugins:
  C compiler: /usr/bin/openmpi/bin/mpicc
  C++ compiler: /usr/bin/openmpi/bin/mpic++
  Fortran compiler: /usr/bin/openmpi/bin/mpif77
  Default backend : CUDA
  CPU backend compiler: /usr/bin/g++ (flags: -w -O3 -g -march=native -ffast-math)
  NVIDIA CUDA backend enabled (flags: -w -O3 -lineinfo --use_fast_math)
  GPU aware MPI support: OFF
  -------------------------------------------
  -- Configuring done (36.8s)
  -- Generating done (0.6s)
  -- Build files have been written to: /home/usr/nekRS/build

  cmake --build ./build --target install -j8
  Please check the summary above carefully and press ENTER to continue or ctrl-c to cancel

Review this output like the installation directory and compilers, and if it looks correct, press **Enter** to compile and install.
On success, you’ll see:

.. code-block:: none

  Hooray! You're all set. The installation is complete.

.. tip::

   Make sure the Default backend matches your hardware: ``CUDA`` (NVIDIA), ``HIP`` (AMD), or ``SYCL`` (Intel).
   If no GPU backend is detected at configure time, NekRS falls back to ``Serial`` (CPU).
   You can switch from a GPU default to Serial at runtime, but not the other way around unless the GPU backend was enabled during build.

Setting the Environment
^^^^^^^^^^^^^^^^^^^^^^^

Assuming the install directory is ``$HOME/.local/nekrs``, set ``NEKRS_HOME`` (required) and optionally prepend its ``bin`` to ``PATH`` (recommended).
For ``bash``:

.. code-block:: bash

   export NEKRS_HOME=$HOME/.local/nekrs
   export PATH=$NEKRS_HOME/bin:$PATH

You may add these lines to ``$HOME/.bashrc`` for persistence.
After the changes, open a new terminal or apply them now:

.. code-block:: bash

   source $HOME/.bashrc

You are now ready to run your *NekRS* cases!


-----------------------------
Running your first simulation
-----------------------------

Several example cases ship with *nekRS* and are installed under ``$NEKRS_HOME/examples``.
To run the ``channel`` example for a simple channel flow case, copy the case folder to scratch location and run it as follows,

.. code-block:: bash

   cd ~
   mkdir scratch
   cd scratch
   cp -r $NEKRS_HOME/examples/channel .
   cd channel
   nrsmpi channel 2

.. note::

   To avoid the conflicts between builds and for cleaner version control, it is recommended to run example case outside both the source directory and install directory.

The final command launches *NekRS* with 2 MPI ranks and print output in the current terminal.
The user may also run *NekRS* in the background as follows:

.. code-block:: bash

   nrsbmpi channel 2

This redirects output to logfiles such as ``channel.log.2`` and creates a ``logfile`` symlink, which users can monitor the progress with:

.. code-block:: bash

   tail -f logfile

.. tip::

   *nekRS* binds 1 GPU to 1 MPI rank. If you saw an OCCA error that fails to create device, you might want to either reduce the number of MPI ranks or fall back to CPU backend via ``nrsbmpi channel 2 --backend serial``.
   The detailed explanation of the argument can be found at :ref:`running`.


-----------------------------
Scripts
-----------------------------

Let’s walk through some useful batch scripts provided by *NekRS*:

- ``nrsmpi <case> #`` runs *NekRS* case in foreground on ``#`` number of processors.
- ``nrsbmpi <case> #`` runs *NekRS* case in background on ``#`` number of processors.
- ``nrsvis <case>`` creates metadata file ``*.nek5000`` required by `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit/>`_ and `ParaView <https://www.paraview.org/>`_. (Usually generated automatically; manual use is rarely needed.)
- ``nrsman env`` lists useful environment variables for *NekRS*
- ``nrsman par`` details all options and settings for the *NekRS* ``.par`` case file

.. _qstart_vis:

-----------------------------
Visualization
-----------------------------

*NekRS* output (``.fld`` or ``0.f%05d``) files can be read by `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit/>`_ or `ParaView <https://www.paraview.org/>`_ via opening the case metadata file ``<case>.nek5000``.
This file is usually created automatically; if missing, generate it with ``nrsvis``.
See :ref:`postproc_checkpoint` for details.
