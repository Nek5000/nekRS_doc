.. _other_resources:

Other Resources
===============

.. _scripts:

Building the Nek5000 Tool Scripts
---------------------------------

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
