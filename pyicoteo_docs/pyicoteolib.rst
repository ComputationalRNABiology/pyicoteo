.. _libdocs:

Pyicoteolib.core
================

Pyicoteolib is the library and the building blocks of the Pyicoteo suite. The pyicoteolib.core library contains the main holders of data in two main objects: ReadCluster and ReadRegion.

ReadCluster
-------------

A ReadCluster object may contain one read or a group of **overlapping** reads, or contigs. It can read both tag like (bed, sam, bam..) and histogram like (wig, bed_pk...) formats. Instances of the ReadCluster object can be added, compared, subtracted to other readCluster objects with standard python syntax.

The ReadCluster object is optimized in order to deal with millions of overlaps, and has been tested with multiple different HTS datasets. The optimization consists in 2 main principles: 

Common python operators
^^^^^^^^^^^^^^^^^^^^^^^^^^

All the following standard operators are supported::

Adding
""""""""

Adding combines the signal of 2 different ReadClusters, with nucleotide precision::

        cluster1 = ReadCluster(read=PK)
        cluster2 = ReadCluster(read=PK)
        cluster1.read_line('chr1 1 45 9:2.00|41:3.00|50:2.00|45:1.00')
        cluster2.read_line('chr1 1 125 9:4.00|41:3.00|30:2.00|45:1.00')
        result = cluster1 + cluster2

        result.write_line()

        chr1    1   145 50:6.00|30:4.00|20:3.00|25:2.00|20:1.00 6.0 .   25  550.0


Subtracting
""""""""""""""

Substracts the signal of 2 different ReadClusters, with nucleotide precision::

        cluster1 = ReadCluster(read=SAM)
        cluster2 = ReadCluster(read=PK)
        cluster1.read_line('SL-XAJ_1_FC305HJAAXX:2:21:872:1402  0   chr1    1   50  36M *   0   0   AAAAGGGGGAATAAAAAGTAACCCAAAACTAACTAT    <<<,7<<<<<7<1:71)<+51<+<5(75()1344+2    PG:Z:FC_305HJAAXX_ln_2.dat')
        cluster2.read_line('chr1 1 125 9:4.00|41:3.00|30:2.00|45:1.00')
        result = cluster2 - cluster1

        result.write_line()

Length
"""""""""

Returns the length of the read cluster::

    c = Cluster(name="chrX", start=1, end=10000)
    len(c)

    10000


Comparison operators (< > == !=)
"""""""""""""""""""""""""""""""""""""

This indicates which read cluster is before another in a chromosome::

    c1 = Cluster(name="chr1", start=100, end=1000)
    c1_copy = Cluster(name="chr1", start=100, end=1000)
    c2 = Cluster(name="chr1", start=50000, end=100000)

    c1 > c2 
    False
    c1 == c1_copy
    True


Lets see some usage examples.

Read a .bed file, print the chromosome and the length of each read
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

    from pyicoteolib.core import ReadCluster, BED

    bed_file = open("/path/to/myfile.bed")

    for line in bed_file:
        rc = ReacCluster(read_as=BED)
        rc.read_line(l)
        print len(rc), rc.area()


Read some .bed lines, cluster them, output a wiggle file
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
::

    cluster =  Cluster(read=BED)
    cluster.read_line('chr1 1 20000 666 hola +')
    cluster.read_line('chr1 1 20000 666 hola +')
    cluster.read_line('chr1 1 20000 666 hola +')
    cluster.read_line('chr1 1001 20000 666 hola +')
    cluster.write_line()


ReadRegion
""""""""""

A ReadRegion object holds a genomic region that may contain ReadClusters

pyicoteolib.utils
------------------

SortedClusterReader
"""""""""""""""""""

SortedCountsReader
"""""""""""""""""""

BigSort
"""""""


BAM reader
------------


Credit
-------

* Developers: Juan González-Vallinas, Ferran Lloret
* Unit and beta Testing: Juan González-Vallinas, Ferran Lloret
* Supervision: Eduardo Eyras


