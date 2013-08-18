
from pyicoteolib.core import ReadRegion, ReadCluster  
from pyicoteolib.utils import SortedFileCountReader, BED

def path_counter(path, region):
	print 

def count_all(paths, experiment_format, region_path, gtf_file, region_magic):
    #first get the region
    if region_path:
    	region_file = open(region_path)

    #load the paths in the readers
    sorted_readers = []
    for path in paths:
        sorted_readers.append(SortedFileCountReader(path, experiment_format))

    #iterate the region objects
    for region_line in region_file:
    	c = ReadCluster(read=experiment_format)
    	c.read_line(region_line)
    	region = ReadRegion(c.name, c.start, c.end)
    	for reader in sorted_readers:
    		print reader.file_path, reader.get_overlaping_counts(region)



