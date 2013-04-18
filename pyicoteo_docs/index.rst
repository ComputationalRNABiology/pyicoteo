.. pyicoteo documentation master file, created by
   sphinx-quickstart on Mon Apr  8 15:27:37 2013.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Pyicoteo
*******************

Pyicoteo* is a suite of tools for the analysis of high-throughput sequencing data. It works with genomic coordinates, it was mainly developed using Solexa/Illumina mapped reads, but it it's core is platform-agnostic. There are currently 6 different tools (5 command-line based, one configuration file based) and a python library for scripting::

	* Pronounced as in Spanish  "picoteo"_ /pɪkɒtɛɒ/: (n) Appetizer-type foods that accompany drinks before or instead of a meal)

If you have any problems or suggestions please join the `Pyicoteo Google Group`_ and ask! 

.. _`Pyicoteo Google Group`: http://groups.google.com/group/pyicos


Getting Started
===============

Download & Install
------------------

**Download Pyicoteo**  `Latest version`_ from our repository.

.. _`Latest version`: https://bitbucket.org/regulatorygenomicsupf/pyicoteo/downloads

You can also download older versions (up to 1.1b) from our `Sourceforge repository`_.

.. _`Sourceforge repository`: http://sourceforge.net/projects/pyicos/ 


The command line tools can be used directly without installation. However, installation is recommended, and necessary if you intend to use the Pyicoteolib. To do so decompress the folder and run the setup.py script with administrator privileges:

    python setup.py install

In order to make it simple for the community, Pyicoteo basic functionality has no dependencies other than Python 2.6 or higher. However, there are 2 optional libraries you could install. 

For plotting capabilities, it is neccesary to install Matplotlib (> 1.0). 

Also, for BAM reading, while we offer a native python implementation, you can ask Pyicoteo to read BAM using samtools with the flag --samtools. 
Pyicoteo is not compatible with Python 3.

Check installation
------------------

To test that the software was installed correctly, start a python console and try importing it by typing::

    python
    >>> import pyicoteolib
    >>> import pyicoteolib.core

The bedpk format
----------------

Some Pyicoteo tools (Pyicos, Pyicaller and Pyicoclip) default experiment and output formats is a derivative of UCSC `Bed format <http://genome.ucsc.edu/FAQ/FAQformat.html#format1>`_ called bedpk. It follows the same starting fields "chromosome/tag start end" but it uses some of the original optional fields to include extra information. It is a cluster oriented format that aims to recollect information of a cluster of reads in a single comprehensive line. 

.. figure:: images/bedpk_format.svg.png 
        :align: left

Column definition
""""""""""""""""""""

1) Chromosome
2) Start coordinate
3) End coordinate
4) Profile: This field summarizes the accumulation of reads per nucleotide of the cluster. The first number is the number of bases covered, while the second will be the number of reads in those bases. See the example above
5) Height: The maximum height of the cluster. In this case, 3.
6) Strand: if ALL clusters reads are positive strand "+", if they are all negative "-". Otherwise "."
7) Summit: The position where the maximum height is found. The binding site is expected to be close to the summit.
8) Area: The area covered by the cluster.
9) p-value: The significance of the cluster calculated by the poisson operation based on peak heights or numbers of reads.

Pyicoteolib
===========

A python library that is the building blocks of all the other tools and a useful tool for python scripting. 

Read more about it at:

.. toctree::
   :maxdepth: 2

	pyicoteolib <pyicoteolib>

Command-line based tools 
========================
.. toctree::
   :maxdepth: 2
  Important considerations <important>
	pyicos <pyicos>
	pyicoller <pyicoller>
	pyicoenrich <pyicoenrich>
	pyicoclip <pyicoclip>
	pyicoregion <pyicoregion>

Protocol files
==============

A configuration file based tool that exposes most functionality of the Pyicoteo suite, making it very useful when trying to combine different tools (for example, Pyicos and Pyicosenrich functionality)

Read more about it at:

.. toctree::
   :maxdepth: 2

	pyicotrocol <pyicotrocol>




