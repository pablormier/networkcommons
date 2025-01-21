############
Installation
############

``NetworkCommons`` requires ``Python`` version >= 3.10 and < 3.13 to run.

-------------------
Requirements
-------------------


NetworkCommons requires `graphviz <https://graphviz.gitlab.io/download/>`_ to visualize networks.

Conda
-------------------

If you are using ``conda``, it can be installed via:

.. code-block:: console

   conda install graphviz

``graphviz`` is also available via the APT package manager in Ubuntu:

Ubuntu
-------------------

.. code-block:: console

   sudo apt-get install -y graphviz graphviz-dev

MacOS
-------------------

In MacOS, it can be installed using `Homebrew <https://brew.sh/>`_. 

.. code-block:: console

   brew install graphviz

In some of our local MacOS tests, we also needed to set the following environment variables:

.. code-block:: console

   export PATH="$(brew --prefix graphviz)/bin:$PATH"'
   export CFLAGS="-I$(brew --prefix graphviz)/include"'
   export LDFLAGS="-L$(brew --prefix graphviz)/lib"'

Please note, if you are using a different architecture or operating system within a subsystem, you may need to make sure that the installed binaries of graphviz were compiled for the correct architecture and that they are visible to the Python interpreter.

-------------------
PIP
-------------------

Currently, NetworkCommons can be installed via pip.

.. code-block:: console

   pip install networkcommons


Additionally, users can install backends for CORNETO via:

.. code-block:: console

   pip install networkcommons[corneto-backends]


To install the development version (with the CORNETO solvers), run:

.. code-block:: console

   pip install 'networkcommons[corneto-backends] @ git+https://github.com/saezlab/networkcommons@dev'


