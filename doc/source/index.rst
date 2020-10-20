rpReader's Documentation
========================

Indices and tables
##################

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

Introduction
############

.. _RetroRules: https://retrorules.org/
.. _RetroPath2.0: https://github.com/Galaxy-SynBioCAD/RetroPath2
.. _RetroPath2.0: https://github.com/Galaxy-SynBioCAD/rp2paths
.. _rpSBML: https://github.com/Galaxy-SynBioCAD/rpBase
.. _rpBase: https://github.com/Galaxy-SynBioCAD/rpBase
.. _rpCache: https://github.com/Galaxy-SynBioCAD/rpCache

Welcome rpReader's documentation. This tool provides a docker that can be accessed using the command line interface for generating rpSBML files from RetroPath2.0_ and rp2paths_, or a TSV file or a string input.

To build the docker you must build a rpBase_ and rpCache_ docker, and then you can use the following command:

.. code-block:: bash

   docker build -t brsynth/rpreader-standalone:v2 .

You can run the docker using the following command to parse the RetroPath2.0_ and rp2paths_ output:

.. code-block:: bash

   python run_rp2.py -rp2paths_compounds test/rp2paths_compounds.csv -rp2_pathways test/rp2_pathways.csv -rp2paths_pathways test/rp2paths_pathways.csv -output test/test_rpReader.tar

If you have a TSV of the pathways (follow the example), you can run the following command:

.. code-block:: bash

   python run_tsv.py -tsv_file /path/to/input_file.tsv -output test/test_rpReader.tar

API
###

.. toctree::
   :maxdepth: 2
   :caption: Contents:

.. currentmodule:: rpTool

.. autoclass:: rpReader
    :show-inheritance:
    :members:
    :inherited-members:

.. currentmodule:: rpToolServe

.. autoclass:: rp2Reader_hdd
    :show-inheritance:
    :members:
    :inherited-members:

.. autoclass:: rp2Reader_mem
    :show-inheritance:
    :members:
    :inherited-members:

.. autoclass:: main_rp2
    :show-inheritance:
    :members:
    :inherited-members:

.. autoclass:: main_tsv
    :show-inheritance:
    :members:
    :inherited-members:

.. autoclass:: main_extrules
    :show-inheritance:
    :members:
    :inherited-members:

.. currentmodule:: run_rp2

.. autoclass:: main
    :show-inheritance:
    :members:
    :inherited-members:

.. currentmodule:: run_tsv

.. autoclass:: main
    :show-inheritance:
    :members:
    :inherited-members:

.. currentmodule:: run_str

.. autoclass:: main
    :show-inheritance:
    :members:
    :inherited-members:
