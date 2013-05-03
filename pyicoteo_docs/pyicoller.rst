Pyicoller
=========

This peak caller is a combination of some of Pyicos commands (extend, normalize, subtract, remove, poisson and filter) for the task of calling peaks from a ChIP-Seq experiment (with narrow peaks). A control file is optional but recommended.


Example::

    pyicoller my_experiment.bed significant_peaks.bedpk -f bed -o --control control.bed --control-format bed --open-control --region regions_to_be_removed.bed --remlabels chrY --correction 0.8 --k-limit 20 --p-value 0.001 -x 130


Credit
------

* Developer: Juan González-Vallinas
* Beta testing: Sonja Althammer, Eneritz Agirre, Nuria Conde Pueyo
* Benchmarking against other peak callers: Sonja Althammer
* Performance benchmarking: Juan González-Vallinas
