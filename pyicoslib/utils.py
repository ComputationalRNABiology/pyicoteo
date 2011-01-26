import math
import os
from tempfile import gettempdir
from heapq import heappop, heappush
from itertools import islice, cycle, chain

from core import Cluster, Region, InvalidLine, InsufficientData, ConversionNotSupported

def add_slash_to_path(path):
    if path[-1] != '/':
        path = '%s/'%path
    return path
    
def poisson(actual, mean):
    '''From StackOverflow: This algorithm is iterative,
        to keep the components from getting too large or small'''
    try:
        p = math.exp(-mean)
        for i in xrange(actual):
            p *= mean
            p /= i+1
        return p
    
    except OverflowError:
        return 0

def pearson(list_one, list_two):
    """
    Accepts paired lists and returns a number between -1 and 1,
    known as Pearson's r, that indicates of how closely correlated
    the two datasets are.
    A score of close to one indicates a high positive correlation.
    That means that X tends to be big when Y is big.
    A score close to negative one indicates a high negative correlation.
    That means X tends to be small when Y is big.
    A score close to zero indicates little correlation between the two
    datasets.
    This script is cobbled together from a variety of sources, linked
    in the sources section below.
    h3. Example usage
    >> import calculate
    >> calculate.pearson([6,5,2], [2,5,6])
    -0.8461538461538467

    h3. A Warning

    Correlation does not equal causation. Just because the two
    datasets are closely related doesn't not mean that one causes
    the other to be the way it is.

    h3. Sources
    http://en.wikipedia.org/wiki/Pearson_product-moment_correlation_coefficient
    http://davidmlane.com/hyperstat/A56626.html
    http://www.cmh.edu/stats/definitions/correlation.htm
    http://www.amazon.com/Programming-Collective-Intelligence-Building-Applications/dp/0596529325
    """
    if len(list_one) != len(list_two):
        raise ValueError('The two lists you provided do not have the name number \
        of entries. Pearson\'s r can only be calculated with paired data.')

    n = len(list_one)

    # Convert all of the data to floats
    list_one = map(float, list_one)
    list_two = map(float, list_two)

    # Add up the total for each
    sum_one = sum(list_one)
    sum_two = sum(list_two)

    # Sum the squares of each
    sum_of_squares_one = sum([pow(i, 2) for i in list_one])
    sum_of_squares_two = sum([pow(i, 2) for i in list_two])

    # Sum up the product of each element multiplied against its pair
    product_sum = sum([item_one * item_two for item_one, item_two in zip(list_one, list_two)])

    # Use the raw materials above to assemble the equation
    pearson_numerator = product_sum - (sum_one * sum_two / n)
    pearson_denominator = math.sqrt((sum_of_squares_one - pow(sum_one,2) / n) * (sum_of_squares_two - pow(sum_two,2) / n))

    # To avoid avoid dividing by zero,
    # catch it early on and drop out
    if pearson_denominator == 0:
        return 0

    # Divide the equation to return the r value
    r = pearson_numerator / pearson_denominator
    return r


def list_available_formats():
    print 'Formats Pyicos can read:'
    for format in READ_FORMATS:
        print format
    print '\nFormats Pyicos can write:'
    for format in WRITE_FORMATS:
        print format
    sys.exit(0)

"""
From Simple Recipes in Python - William Park 1999. 

Returns coefficients to the regression line "y=ax+b" from x[] and y[]. 
""" 

def linreg(X, Y): 
    from math import sqrt   
    if len(X) != len(Y): 
        raise ValueError, 'unequal length' 
    N = len(X) 
    Sx = Sy = Sxx = Syy = Sxy = 0.0 
    for x, y in map(None, X, Y): 
        Sx = Sx + x 
        Sy = Sy + y 
        Sxx = Sxx + x*x 
        Syy = Syy + y*y 
        Sxy = Sxy + x*y 
    det = Sxx * N - Sx * Sx 
    a, b = (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det 

    meanerror = residual = 0.0 
    for x, y in map(None, X, Y): 
        meanerror = meanerror + (y - Sy/N)**2 
        residual = residual + (y - a * x - b)**2 
        RR = 1 - residual/meanerror 
        ss = residual / (N-2) 
        Var_a, Var_b = ss * N / det, ss * Sxx / det 
    return a, b 




class SafeReader:
    def __init__(self, verbose=False):
        self.verbose = verbose
        self.invalid_count = 0
        self.invalid_limit = 2000

    def safe_read_line(self, cluster, line):
        """Reads a line in a file safely catching the error for headers.
        Triggers OperationFailed exception if too many invalid lines are fed to the method"""
        try:
            cluster.read_line(line)
        except InvalidLine:
            if self.invalid_count > self.invalid_limit:
                print
                self.logger.error('Limit of invalid lines: Check the experiment, control, and region file formats, probably the error is in there. Pyicos by default expects bedpk files, except for region files, witch are bed files')
                print
                raise OperationFailed
            else:
                if self.verbose:
                    print "Skipping invalid (%s) line: %s"%(cluster.reader.format, line),
                self.invalid_count += 1
    

class BigSort:
    """
    This class can sort huge files without loading them fully into memory.
    Based on a recipe by Tomasz Bieruta found at http://code.activestate.com/recipes/415581/

    NOTE: This class is becoming a preprocessing module. This is a good thing, I think! But its not
    only a sorting class then. We have to think about renaming it, or extracting functionality from it...
    """
    def __init__(self, file_format=None, read_half_open=False, frag_size=0, id=0, verbose=True):
        self.verbose = verbose
        self.file_format = file_format
        self.frag_size = frag_size
        if self.file_format:
            self.cluster = Cluster(read=self.file_format, write=self.file_format, read_half_open=read_half_open, write_half_open=read_half_open)
        self.id = id
        
    def skipHeaderLines(self, key, experiment_file):
        validLine = False
        count = 0
        while not validLine and count < 40:
            try:
                currentPos = experiment_file.tell()
                line = [experiment_file.readline()]
                line.sort(key=key)
                experiment_file.seek(currentPos)
                validLine = True
            except:
                count += 1

    def filter_chunk(self, chunk):
        filtered_chunk = []
        for line in chunk:
            if self.cluster.reader.quality_filter(line):    
                self.cluster.clear()
                try:           
                    self.cluster.read_line(line)
                    self.cluster.extend(self.frag_size)

                except InvalidLine:
                    print 'Discarding middle invalid line: %s'%line
                                   
                if not self.cluster.is_empty():
                    filtered_chunk.append(self.cluster.write_line())

        return filtered_chunk

    def sort(self, input, output=None, key=None, buffer_size=320000, tempdirs=[], tempFileSize=8000000):
        if key is None:
            key = lambda x : x

        if not tempdirs:
            tempdirs.append(gettempdir())
        input_file = file(input,'rb',tempFileSize)
        self.skipHeaderLines(key, input_file)
        try:
            input_iterator = iter(input_file)
            chunks = []
            for tempdir in cycle(tempdirs):
                current_chunk = list(islice(input_iterator,buffer_size))
                current_chunk = self.filter_chunk(current_chunk) #Now we always filter the chunks, so no empty and invalid lines appear. This used to be optional
                if current_chunk:
                    current_chunk.sort(key=key)
                    output_chunk = file(os.path.join(tempdir,'%06i_%s_%s'%(len(chunks), os.getpid(), self.id)),'w+b',tempFileSize)
                    output_chunk.writelines(current_chunk)
                    output_chunk.flush()
                    output_chunk.seek(0)
                    chunks.append(output_chunk)
                else:
                    break

        except KeyboardInterrupt: #If there is an interruption, delete all temporary files and raise the exception for further processing.
            print 'Removing temporary files...'
            for chunk in chunks:
                try:
                    chunk.close()
                    os.remove(chunk.name)
                except:
                    pass
            raise

        finally:
            input_file.close()
            
        if output is None:
            output = "%s/tempsort%s_%s"%(gettempdir(), os.getpid(), self.id)
        
        output_file = file(output,'wb',tempFileSize)
        
        try:
            output_file.writelines(self.merge(chunks,key))
        finally:
            for chunk in chunks:
                try:
                    chunk.close()
                    os.remove(chunk.name)
                except:
                    pass

        output_file.close()
        return file(output)

    def merge(self, chunks, key=None):
        if self.verbose: print "... Merging chunks..."
        if key is None:
            key = lambda x : x

        values = []
        for index, chunk in enumerate(chunks):
            try:
                iterator = iter(chunk)
                value = iterator.next()
            except StopIteration:
                try:
                    chunk.close()
                    os.remove(chunk.name)
                    chunks.remove(chunk)
                except:
                    pass
            else:
                heappush(values,((key(value),index,value,iterator,chunk)))

        while values:
            k, index, value, iterator, chunk = heappop(values)
            yield value
            try:
                value = iterator.next()
            except StopIteration:
                try:
                    chunk.close()
                    os.remove(chunk.name)
                    chunks.remove(chunk)
                except:
                    pass
            else:
                heappush(values,(key(value),index,value,iterator,chunk))
 

class DualSortedReader:
    """Given two sorted files of tags in a format supported by Pyicos, iterates through them returning them in order"""
    def __init__(self, file_a_path, file_b_path, format, read_half_open=False, verbose=True):
        self.verbose = verbose
        self.file_a = open(file_a_path)
        self.file_b = open(file_b_path)
        self.current_a = Cluster(cached=False, read=format, read_half_open=read_half_open)
        self.current_b = Cluster(cached=False, read=format, read_half_open=read_half_open)
        
    def __iter__(self):
        stop_a = True #indicates if the exception StopIteration is raised by file a (True) or file b (False)
        safe_reader = SafeReader(self.verbose)
        try:
            while 1:
                if not self.current_a:
                    stop_a = True
                    line_a = self.file_a.next()
                    safe_reader.safe_read_line(self.current_a, line_a)
                
                if not self.current_b:
                    stop_a = False
                    line_b = self.file_b.next()
                    safe_reader.safe_read_line(self.current_b, line_b)
                
                if self.current_a < self.current_b:
                    self.current_a.clear()
                    yield line_a
                else:
                    self.current_b.clear()
                    yield line_b
        except StopIteration: #we still need to print the reminder of the sorter file
            if stop_a:
                while self.file_b:
                    yield line_b
                    line_b = self.file_b.next()
            else:
                while self.file_a:
                    yield line_a
                    line_a = self.file_a.next()


class SortedFileClusterReader:
    """
    Holds a cursor and a file path. Given a start and an end, it iterates through the file starting on the cursor position,
    and retrieves the clusters that overlap with the region specified.
    """
    def __init__(self, file_path, experiment_format, read_half_open=False, rounding=True, cached=True):
        self.__dict__.update(locals())
        self.file_iterator = file(file_path)
        self.__initvalues()
        self.safe_reader = SafeReader()
    
    def __initvalues(self):
        self.slow_cursor = 1
        self.cluster_cache = dict() 
        self.invalid_count = 0
        self.invalid_limit = 2000


    def rewind(self):
        """Start again reading the file from the start"""
        self.file_iterator.seek(0)
        self.__initvalues()

    def _read_line_load_cache(self, cursor):
        """Loads the cache if the line read by the cursor is not there yet.
        If the line is empty, it means that the end of file was reached,
        so this function sends a signal for the parent function to halt.
        If the region is stranded, the only tags returned will be the ones of that strand"""
        if cursor not in self.cluster_cache:
            try:
                line = self.file_iterator.next()
            except StopIteration:
                return True
            self.cluster_cache[cursor] = Cluster(read=self.experiment_format, read_half_open=self.read_half_open, rounding=self.rounding, cached=self.cached)
            self.safe_read_line(self.cluster_cache[cursor], line)
        return False

    def get_overlaping_clusters(self, region, overlap=1):
        clusters = []
        if self._read_line_load_cache(self.slow_cursor):
            return clusters
        #advance slow cursor and delete the clusters that are already passed by
        while (self.cluster_cache[self.slow_cursor].chromosome < region.chromosome) or (self.cluster_cache[self.slow_cursor].chromosome == region.chromosome and region.start > self.cluster_cache[self.slow_cursor].end):
            del self.cluster_cache[self.slow_cursor]
            self.slow_cursor+=1
            if self._read_line_load_cache(self.slow_cursor):
                return clusters
        #get intersections
        fast_cursor = self.slow_cursor
        while self.cluster_cache[fast_cursor].start <= region.end and self.cluster_cache[fast_cursor].chromosome == region.chromosome:
            if self.cluster_cache[fast_cursor].overlap(region) >= overlap:
                if not region.strand or region.strand == self.cluster_cache[fast_cursor].strand:
                    clusters.append(self.cluster_cache[fast_cursor].copy_cluster())
            fast_cursor += 1
            if self._read_line_load_cache(fast_cursor):
                return clusters
        return clusters

    def safe_read_line(self, cluster, line):
        self.safe_reader.safe_read_line(cluster, line)

