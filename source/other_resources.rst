.. _other_resources:

Other Resources
===============

.. _linalg:

Linear Algebra Functions (linAlg)
---------------------------------

*TODO*

.. _opsem:

SEM operators (opSEM)
---------------------

*TODO*

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


.. _file_extensions:

File Extensions
---------------

.. csv-table:: Case Files & Extensions
   :header: "Extension","Name","Required?","Notes / Typical usage"
   :widths: 20,30,20,30
   :quote: "
   :delim: ,

   ".sess","NekNek sessions file","Yes (Only for NekNek)","See :ref:`sess_file`"
   ".par","Main parameter file","Yes","See :ref:`par_file`"
   ".udf","User-defined functions","Yes","See :ref:`udf_file`"
   ".oudf","Extra OKL kernels file","No","See :ref:`okl_block`"
   ".upd","Trigger file","No","See :ref:`upd_file`"

   ".usr","Legacy Nek5000 user file","No","See :ref:`usr_file`"
   ".re2","Mesh file","Yes","See :ref:`re2_file`"
   ".co2","Legacy connectivity file","No","See `Nek5000 co2 example <https://nek5000.github.io/NekDoc/tutorials/multi_rans.html?highlight=co2#generating-connectivity-file-co2>`_"

   ".yaml","YAML config file","No","Used for ``nekAscent.hpp`` in-situ visualization"

.. csv-table:: Output Files & Extensions
   :header: "Extension","Name","Notes / Typical usage"
   :widths: 20,40,40
   :quote: "
   :delim: ,

   "\*0.f0000","Nek format checkpoint file","0-based index. See :ref:`qstart_vis`"
   ".nek5000","Nek format checkpoint metadata file", "See: ref: `_qstart_vis`"
   ".bp/","ADIOS2 format checkpoint (folder)",""
   ".log.*","Log files","From ``nrsbmpi``, suffixed by rank count and rollover to ``.log1``"
   ".vtu/.vtk","VTK files","Mainly for particle data"


