pyicoclip
=========

Pyicoclip is an implementation of the modified False Discovery Rate algorithm proposed_ by Yeo et al. to determine which clusters are significant in a list of genomic regions (like genes or transcripts). This method is typically used in CLIP-Seq data that doesn't have a valid control experiment to compare against. 

Theoretically, it could be used for any other kind of experiment that involves short reads and doesn't have a valid control.

A region of interest file is required for the method to be applied, in BED format. 

.. _proposed: http://www.nature.com/nsmb/journal/v16/n2/full/nsmb.1545.html

Example::

    pyicoclip my_experiment.bed my_regions.bed output.pk -f bed

Credit
------

* Developer: Juan González-Vallinas
* Beta Testing: Mireya Plass, Juan González-Vallinas
* Supervision: Eduardo Eyras