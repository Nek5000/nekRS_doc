.. _quickstart:

Quickstart
==========

.. _qstart_before:

-------------------
Before You Begin...
-------------------

It is recommended that the user have basic familiarity with Linux or another Unix-based OS.
The instructions provided in the quickstart guide and the tutorials use basic bash commands and assume the user has this knowledge.
The instructions below are specific to building on your local system (your desktop or laptop). 

.. note::

   Specific instructions for installing on select :term:`HPC` systems can be found `here <https://github.com/Nek5000/nekRS_HPCsupport>`_.

Before you begin installing note the following **requirements** for *nekRS*:

.. mdinclude:: _includes/README.md
   :start-line: 32
   :end-line: 36
   
Most of these should either be available by default in your OS of choice, or can be installed using a common package manager.

.. tabs::

    .. tab:: Debian/Ubuntu

        Debian based systems (such as Ubuntu) use the ``apt`` package manager. GNU
        Compilers for C++ and fortran alongside compatible OpenMPI and CMake
        installs can be acquired with:

        .. code-block:: bash

            sudo apt update
            sudo apt install build-essential libopenmpi-dev cmake

    .. tab:: Mac

        The `Homebrew <https://brew.sh/>`_ package manager is commonly used on Mac
        to provide similar functionality to Linux package managers. GNU Compilers
        for C++ and fortran alongside compatible OpenMPI and CMake installs can be
        acquired with:

        .. code-block:: bash

            brew install gcc open-mpi cmake

**TODO**: Instructions for installing GPU compilers

.. _installing:

----------------
Installing NekRS
----------------

Aquiring the Code
^^^^^^^^^^^^^^^^^

*NekRS* can be acquired by either downloading a release:

.. code-block:: bash

    cd ~
    wget https://github.com/Nek5000/nekRS/archive/refs/tags/v25.0-rc1.tar.gz
    tar -xzvf v25.0-rc1.tar.gz
    mv nekRS-25.0-rc1 nekRS
    cd nekRS

or by cloning the Github repository:

.. code-block:: bash

   cd ~
   git clone https://github.com/Nek5000/nekRS.git
   cd nekRS
   git checkout next

.. note::

   *NekRS* v25 is currently in the ``next`` branch of Github repository.
   Make sure that you checkout the ``next`` branch as shown above.

Building NekRS
^^^^^^^^^^^^^^

From the *NekRS* directory the code can be build using the following general command

.. code-block:: bash

   CC=mpicc CXX=mpic++ FC=mpif77 ./build.sh [-DCMAKE_INSTALL_PREFIX=$HOME/.local/nekrs]

``CC``, ``CXX`` and ``FC`` are environment variables passed to ``cmake`` which specify the C, C++ and fortran compilers, respectively.

.. note::

    On your local system, these environment variables are optional provided you intend to use the system ``MPI`` installed using the ``apt`` or `brew`` package manager as shown :ref:`above <qstart_before>`.
    If the ``MPI`` is installed in a custom location, these must be modified appropriately

    .. code-block:: bash

       CC=/path/to/mpi/bin/mpicc CXX=/path/to/mpi/bin/mpicxx FC=/path/to/mpi/bin/mpif77 ./build.sh

``-DCMAKE_INSTALL_PREFIX`` is a ``cmake`` option which specifies the installation location of *NekRS*.
The default location is ``$HOME/.local/nekrs``.

.. tip::

   Make sure to remove the previous build and installation directory if updating *NekRS*

Running the build script will provide the following summary before compiling *NekRS*

.. code-block:: bash

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

Check the installation directory and compiler info configured by the build system.
Now, to continue to compile and install *NekRS*, just hit enter.

Setting the Environment
^^^^^^^^^^^^^^^^^^^^^^^

Assuming you run ``bash`` and your install directory is ``$HOME/.local/nekrs``, add the following line to your ``$HOME/.bashrc``:

.. code-block:: bash

   export NEKRS_HOME=$HOME/.local/nekrs
   export PATH=$NEKRS_HOME/bin:$PATH

Now, either continue in a new terminal or source the ``$HOME/.bashrc`` file to continue in the current terminal

.. code-block:: bash

   source $HOME/.bashrc

You are now ready to run your *NekRS* cases!

-----------------------------
Running your first simulation
-----------------------------

Several *NekRS* example cases are packaged with the code and provided for the user.
These are located in the ``examples`` folder in the *NekRS* directory.
To run a simple channel flow case included in this folder, copy the example case to another location and run it as follows,

.. code-block:: bash

   cd ~
   mkdir scratch
   cd scratch
   cp -r ~/nekRS/examples/channel .
   cd channel
   nrsmpi channel 2

.. note::

   It is recommended that the example case should be copied to another directory as shown above before running.
   Avoid running in the ``nekRS/example`` folder.

The final command will launch *NekRS* and the output will be shown in the current terminal window.
The user may also run *NekRS* in the background as follows,

.. code-block:: bash

   nrsbmpi channel 2

The above command will print output in ``logfile`` which the user can monitor as the case in running in the background using

.. code-block:: bash

   tail -f logfile

-----------------------------
Scripts
-----------------------------

Let’s walk through some useful batch scripts provided by *NekRS*:

- ``nrsmpi <case> #`` runs *NekRS* case in foreground on ``#`` number of processors.
- ``nrsbmpi <case> #`` runs *NekRS* case in background on ``#`` number of processors.
- ``nrsvis <case>`` creates metadata file required by `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit/>`_ and `ParaView <https://www.paraview.org/>`_.
- ``nrsman env`` lists useful environment variables for *NekRS*
- ``nrsman par`` details all options and settings for the *NekRS* ``.par`` case file

-----------------------------
Visualization
-----------------------------

*NekRS* output (``.fld`` or ``0.f%05d``) files can be read by `VisIt <https://wci.llnl.gov/simulation/computer-codes/visit/>`_ or `ParaView <https://www.paraview.org/>`_.
This requires using ``nrsvis <case>`` to generate a metadata file.
