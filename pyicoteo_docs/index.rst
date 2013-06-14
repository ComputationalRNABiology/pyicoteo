.. pyicoteo documentation master file, created by
   sphinx-quickstart on Mon Apr  8 15:27:37 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Pyicoteo
*******************

Pyicoteo* is a suite of tools for the analysis of high-throughput sequencing data. It works with genomic coordinates, it was mainly developed using Solexa/Illumina mapped reads, but it it's core is platform-agnostic. There are currently 6 different tools for the analysis of HTS data: 5 command-line based, one configuration file based and a python library for scripting::

	* Pronounced as in Spanish  "picoteo"_ /pɪkɒtɛɒ/: 
          (n) Appetizer-type foods that accompany drinks before or instead of a meal.


Please start by reading the :ref:`intro` document. 

Protocol files
==============

A configuration file based tool that exposes most functionality of the Pyicoteo suite, making it very useful when trying to combine different tools (for example, Pyicos and Pyicoenrich functionality)

Read more about it at :ref:`protocoldocs`


Full Table of Contents
======================

.. toctree::
  :maxdepth: 2

  intro
  pyicos
  pyicoller
  pyicoenrich
  pyicoclip
  pyicoregion
  pyicotrocol
  pyicoteolib

