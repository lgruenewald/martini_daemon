Martini Daemon documentation
============================

Martini Daemon is a tool facilitating template based chemical reactions in MD simulations with the
`Martini force field`_ and the `OpenMM`_ MD engine.

.. _Martini force field: https://cgmartini.nl/
.. _OpenMM: https://openmm.org/

This is achieved by combining multiple components in one python package:

* a friendly Python API for running MD simulations with reactions,
* a parser for GROMACS ``.top`` files targeting OpenMM,
* a thin wrapper on top of OpenMM's API facilitating bond addition and removal,
* a graph matching system to find reactants,
* a detection/modification algorithm to execute reaction templates.

.. Note::

    See the README of the git repository for installation instructions.

    After installing Martini Daemon, new users are advised to read the :doc:`/tutorial` first.

.. toctree::
   :maxdepth: 2
   :caption: Contents:

   tutorial
   user_guide
   faq
   reference
   extending

