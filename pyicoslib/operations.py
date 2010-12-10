"""
Pyicos is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

import os.path
import sys
import os
import shutil
import logging
import math
import random
from collections import defaultdict
from heapq import heappop, heappush
from itertools import islice, cycle
from tempfile import gettempdir
from datetime import datetime
from core import Cluster, Region, InvalidLine, InsufficientData, ConversionNotSupported
from defaults import *
from tempfile import gettempdir

Normalize = 'normalize'
Extend = 'extend'
Subtract = 'subtract'
Split = 'split'
Trim = 'trim'
Cut = 'filter'
Poisson = 'poisson'
Convert = 'convert'
NoWrite = 'nowrite'
DiscardArtifacts = 'discard'
RemoveRegion = 'remove'
RemoveDuplicates = 'remove_duplicates'
ModFDR = 'modfdr'
StrandCorrelation = 'strand_correlation'
Enrichment = 'enrichment'


class OperationFailed(Exception):
    pass

class Utils:
    @staticmethod
    def add_slash_to_path(path):
        if path[-1] != '/':
            path = '%s/'%path
        return path
    

    @staticmethod
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

    @staticmethod
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

    @staticmethod
    def list_available_formats():
        print 'Formats Pyicos can read:'
        for format in READ_FORMATS:
            print format
        print '\nFormats Pyicos can write:'
        for format in WRITE_FORMATS:
            print format
        sys.exit(0)

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

        def sort(self, input,output=None,key=None,buffer_size=320000, tempdirs=[], tempFileSize=8000000):
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

        def merge(self, chunks,key=None):
            if self.verbose: print "Merging chunks..."
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
                    
                    
class SortedFileClusterReader:
    """Holds a cursor and a file path. Given a start and an end, it iterates through the file starting on the cursor position,
    and retrieves the clusters that overlap with the region specified.
    """
    def __init__(self, file_path, experiment_format, read_half_open=False, rounding=True, cached=True):
        self.__dict__.update(locals())
        self.slow_cursor = 1
        self.cluster_cache = dict() 
        self.invalid_count = 0
        self.invalid_limit = 2000
        self.file_iterator = file(file_path)

    def _read_line_load_cache(self, cursor):
        """Loads the cache if the line read by the cursor is not there yet.
        If the line is empty, it means that the end of file was reached,
        so this function sends a signal for the parent function to halt """
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
                clusters.append(self.cluster_cache[fast_cursor].copy_cluster())
            fast_cursor += 1
            if self._read_line_load_cache(fast_cursor):
                return clusters
        return clusters

    #TODO read_safe_line dentro de Cluster, o Utils...
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
                print "Skipping invalid (%s) line: %s"%(cluster.reader.format, line),
                self.invalid_count += 1



class Turbomix:
    """
    This class is the pipeline that makes possible the different combination of operations. 
    It has different switches that are activated by the list 'self.operations'.
    """
    logger = logging.getLogger("log/pyicos.log")
    logger.setLevel(logging.WARNING)
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)
    formatter = logging.Formatter("%(levelname)s - %(message)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    invalid_count = 0
    invalid_limit = 2000
    control_path = None
    
    def __init__(self, experiment_path, output_path, experiment_format=BED, output_format=PK, label=LABEL, 
                 open_experiment=OPEN_EXPERIMENT, open_output=OPEN_OUTPUT, debug = DEBUG, rounding=ROUNDING, tag_length = TAG_LENGTH, discarded_chromosomes = REMLABELS,
                 control_path = CONTROL, control_format = PK, open_control = OPEN_CONTROL, region_path = REGION, region_format = PK, 
                 open_region=OPEN_REGION, span = SPAN, frag_size = FRAG_SIZE, p_value = P_VALUE, height_limit = HEIGHT_LIMIT, correction_factor = CORRECTION,
                 trim_percentage=TRIM_PROPORTION, no_sort=NO_SORT, duplicates=DUPLICATES, threshold=THRESHOLD, trim_absolute=TRIM_ABSOLUTE, max_delta=MAX_DELTA, 
                 min_delta=MIN_DELTA, correlation_height=HEIGHT_FILTER, delta_step=DELTA_STEP, verbose=VERBOSE, species=SPECIES, cached=CACHED, split_proportion=SPLIT_PROPORTION, 
                 split_absolute=SPLIT_ABSOLUTE, repeats=REPEATS, masker_file=MASKER_FILE, max_correlations=MAX_CORRELATIONS, keep_temp=KEEP_TEMP, experiment_b_path=EXPERIMENT,
                 replica_a_path=EXPERIMENT, replica_b_path=EXPERIMENT, poisson_test=POISSONTEST, **kwargs):
        self.__dict__.update(locals())
        self.is_sorted = False
        self.temp_experiment = False #Indicates if temp files where created for the experiment
        self.temp_control = False #Indicates if temporary files where created for the control
        self.temp_replica_a = False 
        self.temp_replica_b = False 
        self.operations = []
        self.previous_chr = None
        self.open_region = False
        
        if not self.discarded_chromosomes:
            self.discarded_chromosomes = []
        self.region_cluster = Cluster(read=region_format, read_half_open = open_region)
        #Reusable cluster objects 
        try:
            self.cluster = Cluster(read=self.experiment_format, write=self.output_format, rounding=self.rounding, read_half_open = self.open_experiment, write_half_open = self.open_output, tag_length=self.tag_length, span = self.span, verbose=self.verbose, cached=cached)
            self.cluster_aux = Cluster(read=self.experiment_format, write=self.output_format, rounding=self.rounding, read_half_open = self.open_experiment, write_half_open = self.open_output, tag_length=self.tag_length, span = self.span, verbose=self.verbose, cached=cached)
            self.cluster_aux2 = Cluster(read=self.experiment_format, write=self.output_format, rounding=self.rounding, read_half_open = self.open_experiment, write_half_open = self.open_output, tag_length=self.tag_length, span = self.span, verbose=self.verbose, cached=cached)
            self.control_cluster = Cluster(read=experiment_format, read_half_open = open_experiment, verbose=self.verbose, cached=cached)

        except ConversionNotSupported:
            print '\nThe reading "%s" and writing "%s" is not supported. \n\n'%(self.experiment_format, self.output_format)
            Utils.list_available_formats()
        #duplicates flag
        self.previous_start = 0
        self.previous_end = 0
        self.previous_chr = ''
        self.duplicates_found = 0
        #poisson stuff
        self.first_chr = True
        self._init_poisson()
        self.poisson_results = {'length': defaultdict(), 'height': defaultdict(), 'numreads': defaultdict()}
        self.maxheight_to_pvalue = {}
        self.numreads_to_pvalue = {}
        #Operation flags
        self.do_poisson = False
        self.do_subtract = False
        self.do_heuremove = False
        self.do_split = False
        self.do_trim = False
        self.do_cut = False
        self.do_extend = False
        self.do_dupremove = False
        self.sorted_by_pyicos = False
        self.sorted_region_path = ''

    def i_cant_do(self):
        """Quits and exists if exits if the combination of non possible operations"""
        if Cut in self.operations and Poisson not in self.operations and not self.threshold:
            print "Can't do Cut without Poisson or a fixed threshold\n"
            sys.exit(0)
        elif (Subtract in self.operations or Split in self.operations or Poisson in self.operations) and (self.output_format not in CLUSTER_FORMATS):
            print 'Cant get the output as tag format (eland, bed) for Subtract, Split, Poisson filtering please use a clustered format %s\n'%CLUSTER_FORMATS
            sys.exit(0)
        elif Extend in self.operations and self.experiment_format in CLUSTER_FORMATS:
            print "Can't extend if the experiment is a clustered format ",
            print CLUSTER_FORMATS
            sys.exit(0)
        elif StrandCorrelation in self.operations and self.experiment_format in CLUSTER_FORMATS:
            print "Can't perform strand correlation operation if the experiment is a clustered format ",
            print CLUSTER_FORMATS
            sys.exit(0)
    
        
    def _add_slash_to_path(self, path):
        return Utils.add_slash_to_path(path)

    def read_and_preprocess(self, cluster, line):
        self.safe_read_line(cluster, line)
        if not cluster.is_empty():
            if self.do_heuremove:
                cluster = self.remove_regions(cluster)
            if self.do_extend:
                cluster.extend(self.frag_size)


    def success_message(self, output_path):
        self.result_log.write('\nSummary of operations\n')
        self.result_log.write('---------------------\n')
        if self.verbose:
            print 'Success!\n'
            print 'Summary of operations:'
        
            if self.experiment_format != self.output_format and not NoWrite in self.operations:
                print 'Convert from %s to %s.'%(self.experiment_format, self.output_format)

        if StrandCorrelation in self.operations:
            self.result_log.write('Strand correlation: %.0f pairs analyzed, estimated a %s nucleotides extension value'%(self.analyzed_pairs, self.frag_size))
            if self.verbose: print 'Strand correlation'

        if Extend in self.operations:
            self.result_log.write('Extend to %s\n'%self.frag_size)
            if self.verbose: print 'Extend to %s'%self.frag_size

        if self.do_subtract:
            self.result_log.write('Subtract\n')
            if self.verbose: print 'Subtract'

        if self.discarded_chromosomes:
            self.result_log.write('Discard tags: %s\n'%self.discarded_chromosomes)
            if self.verbose: print 'Discard tags: %s'%self.discarded_chromosomes

        if self.do_heuremove:
            self.result_log.write('Heuristic Remove from %s\n'%self.region_path)
            if self.verbose: print 'Heuristic Remove from %s'%self.region_path

        if self.do_split:
            self.result_log.write('Split\n')
            if self.verbose: print 'Split'

        if DiscardArtifacts in self.operations:
            self.result_log.write('Discard Artifacts\n')
            if self.verbose: print 'Discard Artifacts'

        if self.do_poisson:
            self.result_log.write('Poisson\n')
            if self.verbose: print 'Poisson'

        if self.do_cut:
            self.result_log.write('Filter\n')
            if self.verbose: print 'Filter'

        if self.do_dupremove:
            self.result_log.write('Removed %s duplicates\n'%self.duplicates_found)
            if self.verbose: print 'Removed %s duplicates, allowing only up to %s'%(self.duplicates_found, self.duplicates)

        self.result_log.write('Date finished: %s'%datetime.now())        
        if not NoWrite in self.operations and self.verbose:
            print 'Output at: %s'%(output_path)

    def start_operation_message(self):
        if self.verbose:
            if self.experiment_format != self.output_format and not NoWrite in self.operations:
                print 'Converting file %s from %s to %s...'%(self.current_experiment_path, self.experiment_format, self.output_format)
            else:
                print 'Reading file %s as a %s file...'%(self.current_experiment_path, self.experiment_format)

            if self.current_control_path:
                print 'Control file:%s'%self.current_control_path
                if Normalize in self.operations:
                    print 'The file %s will be normalized to match %s'%(self.current_experiment_path, self.current_control_path)
                if Subtract in self.operations:
                    print 'The file %s will be subtracted from %s'%(self.current_control_path, self.current_experiment_path)

    def get_normalize_factor(self, experiment, control):
        ret = self.numcells(control, self.control_format)/self.numcells(experiment, self.experiment_format)
        self.result_log.write('Normalization factor: %s\n'%(ret))
        print 'Normalization factor: %s\n'%(ret)
        return ret
    
    def numcells(self, file_path, file_format):
        """Returns the total number of cells in a file"""
        num = 0.
        cluster = Cluster(read=file_format)
        for line in file(file_path, 'rb'):
            self.safe_read_line(cluster, line)
            for level in cluster: 
                num += level[0]*level[1]
            cluster.clear()

        return num  
    
    def safe_read_line(self, cluster, line):
        """Reads a line in a file safely catching the error for headers. 
        Triggers OperationFailed exception if too many invalid lines are fed to the method"""
        try:
            cluster.read_line(line)
        except InvalidLine:
            if self.invalid_count > self.invalid_limit:
                print 
                self.logger.error('Limit of invalid lines: Incorrect file format? Check the experiment, control, and region formats, probably the error is in there. Pyicos by default expects bedpk files, or regular bed files for the region files.')
                print
                Utils.list_available_formats()
                raise OperationFailed
            else:
                #if self.verbose: print "Skipping invalid (%s) line: %s"%(cluster.reader.format, line),
                self.invalid_count += 1

    def run(self):
        self.do_subtract = (Subtract in self.operations and self.control_path is not None)
        self.do_normalize = (Normalize in self.operations and self.control_path is not None)
        self.do_heuremove = (RemoveRegion in self.operations and self.region_path)
        self.do_poisson = Poisson in self.operations
        self.do_split = Split in self.operations
        self.do_trim = Trim in self.operations
        self.do_cut = Cut in self.operations
        self.do_extend = Extend in self.operations and not self.sorted_by_pyicos #If the experiment was sorted by Pyicos, it was already extended before, so dont do it again
        self.do_discard = DiscardArtifacts in self.operations
        self.do_dupremove = RemoveDuplicates in self.operations
        if self.verbose:
            print '\n\nPyicos running...'
        if self.control_path:
            self.process_all_files_paired(self.experiment_path, self.control_path)
        elif self.experiment_b_path:
            self.process_all_files_paired(self.experiment_path, self.experiment_b_path)
        else:
            self.process_all_files_recursive(self.experiment_path, self.output_path)

    def process_all_files_paired(self, path_a, path_b):
        if os.path.isdir(path_a):
            self.logger.warning('Operating with directories will only give appropiate results if the files and the control are paired in alphabetical order.')
            controls = os.listdir(path_b)
            controls.sort()
            experiments = os.listdir(path_a)
            experiments.sort()
            for i in xrange(0, len(experiments)):
                try:
                    operate(experiments[i], control[i], '%s%s'%(self.output_path, experiments[i]))
                except IndexError:
                    pass
        else:
            output = self.output_path
            if os.path.isdir(self.output_path):
                output = self._add_slash_to_path(self.output_path)
                output = '%s%s_minus_%s'%(self.output_path, os.path.basename(path_a), os.path.basename(path_b))
            
            self.operate(path_a, path_b, output)
    
    def process_all_files_recursive(self, dirorfile, output, firstIteration=True):
        paths = []
        """Goes through all the directories and creates files recursively. Returns the paths of the resulting files"""
        if os.path.isdir(dirorfile):
            dirorfile = self._add_slash_to_path(dirorfile)
            if not firstIteration:
                output = '%s/%s/'%(os.path.abspath(output), dirorfile.split('/')[-2])
                if not os.path.exists(output):
                    os.makedirs(output)  
            
            for filename in os.listdir(dirorfile):
                self.process_all_files_recursive('%s%s'%(dirorfile,filename), output)
                
        elif os.path.isfile(os.path.abspath(dirorfile)):
            if os.path.isdir(output):
                output = '%s%s.%s'%(self._add_slash_to_path(output), (os.path.basename(dirorfile)).split('.')[0], self.output_format)
            try:
                self.operate(experiment_path=os.path.abspath(dirorfile), output_path=os.path.abspath(output))
                paths.append(output)
            except OperationFailed:
                self.logger.warning('%s couldnt be read. Skipping to next file.'%os.path.abspath(dirorfile))
                os.remove(os.path.abspath(output))
            except StopIteration:
                self.logger.warning('%s End of file reached.'%os.path.abspath(dirorfile))
        else:
            self.logger.error('Input "%s" doesnt exist?'%os.path.abspath(dirorfile))

        return paths

    def normalize(self):
        if self.control_path:
            if self.verbose: print 'Calculating normalization factor...'
            self.normalize_factor = self.get_normalize_factor(self.current_experiment_path, self.current_control_path)
            self.cluster.normalize_factor = self.normalize_factor
            self.cluster_aux.normalize_factor = self.normalize_factor


    def _manage_temp_file(self, path):
        """A temporary file that is no longer needed is given, and depending on the value of self.keep_temp its removed or kept"""
        if self.keep_temp:
            print "Temp file kept at: %s"%path
        else:
            os.remove(path)
            if self.verbose:
                print 'Temporary file %s removed'%path

    def _remove_duplicates_file(self, file_path, temp_name, remove_temp):
        new_file_path = "%s/tempfilter%s_%s"%(gettempdir(), os.getpid(), temp_name)
        new_file = file(new_file_path, 'w')
        print "Filtering %s file..."%temp_name
        previous_line = ''
        equal_lines = 0
        self.cluster.clear()
        self.cluster_aux.clear()
        for line in file(self.current_experiment_path):
            self.cluster.clear()
            self.cluster.read_line(line)
            if self.cluster.start == self.cluster_aux.start and self.cluster.end == self.cluster.end and self.cluster.chromosome == self.cluster.chromosome:
                equal_lines+=1
            else:
                equal_lines=0
            self.cluster_aux.clear()
            self.cluster_aux.read_line(line)
           
            if self.duplicates >= equal_lines:
                new_file.write(line)
            else:
                self.duplicates_found += 1

        new_file.flush()
        new_file.close()
        if remove_temp:
            self._manage_temp_file(file_path) 
    
        self.cluster.clear()
        self.cluster_aux.clear()   
        return new_file_path     
    
    def filter_files(self):
        """Filter the files removing the duplicates, calculates the enrichment based on the reads."""
        #TODO Normalize should go here too
        if self.do_dupremove:
            self.current_experiment_path = self._remove_duplicates_file(self.current_experiment_path, "experiment", self.temp_experiment)
            if self.current_control_path:
                self.current_control_path = self._remove_duplicates_file(self.current_control_path, "control", self.temp_control)


    def get_lambda_func(self, format):
        if self.output_format == SPK:
            if format == ELAND:
                return lambda x:(x.split()[6],x.split()[8],int(x.split()[7]),len(x.split()[1]))
            else:
                return lambda x:(x.split()[0],x.split()[5],int(x.split()[1]),int(x.split()[2]))
        else:
            if format == ELAND:
                return lambda x:(x.split()[6], int(x.split()[7]), len(x.split()[1]))
            elif format == SAM:
                return lambda x:(x.split()[2], int(x.split()[3]), len(x.split()[9]))
            else:
                return lambda x:(x.split()[0],int(x.split()[1]),int(x.split()[2]))


    def decide_sort(self, experiment_path, control_path=None):
        """Decide if the files need to be sorted or not."""
        #TODO refractor this, copy pasted code (not as easy as it seems)

        if (not self.experiment_format in CLUSTER_FORMATS and self.output_format in CLUSTER_FORMATS) or self.do_subtract or self.do_heuremove or self.do_dupremove or ModFDR in self.operations or Enrichment in self.operations:
            self.experiment_preprocessor = Utils.BigSort(self.experiment_format, self.open_experiment, self.frag_size, 'experiment', verbose=self.verbose)
            self.experiment_b_preprocessor = Utils.BigSort(self.experiment_format, self.open_experiment, self.frag_size, 'experiment_b', verbose=self.verbose)
            self.replica_a_preprocessor = Utils.BigSort(self.experiment_format, self.open_experiment, self.frag_size, 'replica_a', verbose=self.verbose)
            self.replica_b_preprocessor = Utils.BigSort(self.experiment_format, self.open_control, self.frag_size, 'replica_b', verbose=self.verbose)
            self.control_preprocessor = Utils.BigSort(self.control_format, self.open_control, self.frag_size, 'control', verbose=self.verbose)
            self.region_preprocessor = Utils.BigSort(self.region_format, self.open_region, None, 'region', verbose=self.verbose)
            if self.no_sort:
                if self.verbose:
                    print 'Input sort skipped'
                self.current_experiment_path = experiment_path
            else:
                if self.verbose: print 'Sorting experiment file...'
                self.is_sorted = True
                sorted_experiment_file = self.experiment_preprocessor.sort(experiment_path, None, self.get_lambda_func(self.experiment_format))
                self.current_experiment_path = sorted_experiment_file.name
                self.temp_experiment = True

            if self.do_subtract:
                if self.no_sort:
                    if self.verbose: print 'Control sort skipped'
                    self.current_control_path = control_path
                else:
                    if self.verbose: print 'Sorting control file...'
                    sorted_control_file = self.control_preprocessor.sort(control_path, None, self.get_lambda_func(self.control_format))
                    self.current_control_path = sorted_control_file.name
                    self.temp_control = True
            
            if Enrichment in self.operations:  
                if self.no_sort:
                    if self.verbose: print 'Experiment_b sort skipped'
                    self.current_control_path = control_path
                else:
                    if self.verbose: print 'Sorting experiment_b file...'
                    sorted_control_file = self.experiment_b_preprocessor.sort(control_path, None, self.get_lambda_func(self.experiment_format))
                    self.current_control_path = sorted_control_file.name
                    self.temp_control = True
            
            if self.replica_a_path:
                if self.no_sort:
                    if self.verbose: print 'replica_a sort skipped'
                    self.current_replica_a_path = self.replica_a_path
                else:
                    if self.verbose: print 'Sorting replica_a file...'
                    sorted_replica_a_file = self.replica_a_preprocessor.sort(self.replica_a_path, None, self.get_lambda_func(self.experiment_format))
                    self.current_replica_a_path = sorted_replica_a_file.name
                    self.temp_replica_a = True

            if self.replica_b_path:
                if self.no_sort:
                    if self.verbose: print 'replica_b sort skipped'
                    self.current_replica_b_path = self.replica_b_path
                else:
                    if self.verbose: print 'Sorting replica_b file...'
                    sorted_replica_b_file = self.replica_b_preprocessor.sort(self.replica_b_path, None, self.get_lambda_func(self.experiment_format))
                    self.current_replica_b_path = sorted_replica_b_file.name
                    self.temp_replica_b = True


               
        if self.region_path:
            print "Sorting region file..."
            self.sorted_region_file = self.region_preprocessor.sort(self.region_path, None, self.get_lambda_func(BED))
            self.sorted_region_path = self.sorted_region_file.name

    def operate(self, experiment_path, control_path=None, output_path=None):
        """Operate expects single paths, not directories. Its called from run() several times if the experiment for picos is a directory"""
        try:
            self.i_cant_do()
            #per operation variables
            self.previous_chr = None
            self.current_experiment_path = experiment_path
            self.current_control_path = control_path
            self.current_output_path = output_path
            self.cluster.clear()
            self.cluster_aux.clear()
            self.result_log = file('%s/pyicos_report_%s.txt'%(os.path.dirname(os.path.abspath(self.current_output_path)), os.path.basename(self.current_experiment_path)), 'wb')
            self.result_log.write('Pyicos analysis report\n')
            self.result_log.write('----------------------\n\n')
            self.result_log.write('Date run: %s\n'%datetime.now())
            self.result_log.write('Experiment file: %s\n'%self.current_experiment_path)
            if self.current_control_path:
                self.result_log.write('Control file: %s\n'%self.current_control_path)
            self.result_log.write('\n\n')
            if self.output_format == WIG and self.open_output == False:
                self.logger.warning('You are creating a closed wig file. This will not be visible in the UCSC genome browser')

            self.start_operation_message()
            self.decide_sort(experiment_path, control_path)

            self.estimate_frag_size = self.do_poisson and not self.frag_size

            if self.current_control_path:
                self.control_reader = SortedFileClusterReader(self.current_control_path, self.control_format, rounding=self.rounding, cached=self.cached)
            
            if self.region_path:
                self.region_reader = SortedFileClusterReader(self.sorted_region_path, BED, rounding=self.rounding, cached=self.cached)

            if StrandCorrelation in self.operations:
                self.strand_correlation()
            
            if self.do_dupremove:
                self.filter_files()

            if Normalize in self.operations:
                self.normalize()

            if self.do_cut: #if we cut, we will round later
                self.cluster.rounding = False
                self.cluster_aux.rounding = False

            if ((not NoWrite in self.operations) or (NoWrite in self.operations and Poisson in self.operations)) and not Enrichment in self.operations: 
                self.process_file()
            
            if self.do_poisson: #extract info from the last chromosome and print the thresholds
                self.poisson_analysis(self.previous_chr)
                print '\nCluster threadsholds for p-value %s:'%self.p_value
                self.result_log.write('\nCluster threadsholds for p-value %s:\n'%self.p_value)
                for chromosome, k in self.poisson_results[self.poisson_test].items():
                    print '%s: %s'%(chromosome, k)
                    self.result_log.write('%s: %s\n'%(chromosome, k))


            #Mutually exclusive final operations
            if Enrichment in self.operations:
                self.enrichment()

            elif self.do_cut: 
                self.cut()

            elif ModFDR in self.operations:
                self.modfdr()

            self.success_message(output_path)

            

        finally: #Finally try deleting all temporary files, quit silently if they dont exist
            try:
                if self.temp_experiment:
                    self._manage_temp_file(self.current_experiment_path)

            except AttributeError, OSError:
                pass
            try:
                if self.temp_control:
                    self._manage_temp_file(self.current_control_path)

            except AttributeError, OSError:
                pass
 
            try:
                if self.sorted_region_path:
                    self._manage_temp_file(self.sorted_region_file.name)

            except AttributeError, OSError:
                pass


    def _to_read_conversion(self, experiment, output):
        for line in experiment:
            try:
                self.read_and_preprocess(self.cluster, line)
                if not self.cluster.is_empty():
                    self.process_cluster(self.cluster, output)
                self.cluster.clear()
            except InsufficientData:
                self.logger.warning('For this type of conversion (%s to %s) you need to specify the tag length with the --tag-length flag'%(self.experiment_format, self.output_format))
                sys.exit(0)


    def _to_cluster_conversion(self, experiment, output):
        while self.cluster.is_empty():
            self.cluster_aux2.clear()
            self.read_and_preprocess(self.cluster_aux2, experiment.next())
            if not self.cluster_aux2.is_empty():
                self.cluster = self.cluster_aux2.copy_cluster()
            
        for line in experiment:
            self.cluster_aux.clear()
            self.read_and_preprocess(self.cluster_aux, line)
            if not self.cluster_aux.is_empty():
                if self.cluster_aux2.intersects(self.cluster_aux) or self.cluster_aux.is_contiguous(self.cluster_aux2) or self.cluster_aux2.is_contiguous(self.cluster_aux):
                    self.cluster += self.cluster_aux 
                    self.cluster_aux2.clear()
                    self.cluster_aux2 = self.cluster_aux.copy_cluster()
                else:
                    if not self.cluster.is_empty():
                        self.process_cluster(self.cluster, output)
                    self.cluster.clear()
                    self.cluster_aux2.clear()
                    self.read_and_preprocess(self.cluster_aux2, line)
                    if not self.cluster_aux2.is_empty():
                        self.cluster = self.cluster_aux2.copy_cluster()  

        if not self.cluster.is_empty():
            self.process_cluster(self.cluster, output)
            self.cluster.clear()


    def process_file(self):
        self.cluster.clear()
        self.cluster_aux.clear()
        if NoWrite in self.operations:
            output = None
        else:
            output = file(self.current_output_path, 'wb')
        experiment = file(self.current_experiment_path, 'rb')

        if self.output_format == WIG or self.output_format == VARIABLE_WIG:
            output.write('track type=wiggle_0\tname="%s"\tvisibility=full\n'%self.label)

        if self.output_format in CLUSTER_FORMATS and not ModFDR in self.operations and not Enrichment in self.operations:
            if not self.experiment_format in CLUSTER_FORMATS and self.verbose:
                print 'Clustering reads...'
            self._to_cluster_conversion(experiment, output)
        else:
            self._to_read_conversion(experiment, output)

        if not NoWrite in self.operations:
            output.flush()
            output.close()

    def subtract(self, cluster):
        region = Region(cluster.start, cluster.end, cluster.name, cluster.chromosome)
        over = self.control_reader.get_overlaping_clusters(region, overlap=0.0000001)
        for tag in over:
            if self.do_extend:
                tag.extend(self.frag_size)
        region.add_tags(over, True)
        meta = region.get_metacluster()
        cluster -= meta
        return cluster

    def remove_regions(self, cluster):
        region = Region(cluster.start, cluster.end, cluster.name, cluster.chromosome)
        if self.region_reader.get_overlaping_clusters(region, overlap=0.5):
            cluster.clear() 
        return cluster

    def remove_duplicates(self, cluster):
        """Removes the duplicates found in the experiment file"""
        if cluster.start == self.previous_start and self.previous_end == cluster.end and self.previous_chr == cluster.chromosome and not cluster.is_empty():
            self.identical_reads+=1
            if self.identical_reads > self.duplicates:
                self.duplicates_found += 1
                cluster.clear()
        else:
            self.previous_start = cluster.start
            self.previous_end = cluster.end
            self.previous_chr = cluster.chromosome
            self.identical_reads = 0

        return cluster

    def _correct_bias(self, p_value):
        if p_value < 0:
            return 0
        else:
            return p_value

    def _init_poisson(self):
        self.total_bp_with_reads = 0.
        self.total_clusters = 0.
        self.total_reads = 0.
        self.acum_height = 0.
        self.absolute_max_height = 0
        self.absolute_max_numreads = 0
        self.chr_length = 0

    def poisson_analysis(self, chromosome=''):
        """
        We do 3 different poisson statistical tests per chromosome for each experiment file:
        
        Nucleotide analysis:
        This analysis takes the nucleotide as the unit to analize. We give a p-value for each "height"
        of read per nucleotides using an accumulated poisson. With this test we can infer more accurately 
        what nucleotides in the cluster are part of the DNA binding site.
 
        Cluster analysis:
        This analysis takes as a basic unit the "cluster" profile and performs a poisson taking into account the height
        of the profile. This will help us to better know witch clusters are statistically significant and which are product of chromatine noise

        We do this by calculating the acumulated maximum height for all existing clusters and dividing it by the number of clusters. This gives us the average height of a cluster.
        Given a mean and a   height k the poisson function gives us the probability p of one cluster having a height k by chance. With this if we want, for example,
        to know what is the probability of getting a cluster higher than k = 7, we accumulate the p-values for poisson(k"0..6", mean).
        
        Number of reads analysis:
        We analize the number of reads of the cluster. Number of reads = sum(xi *yi ) / read_length

        """
        if os.path.isdir(self.output_path):
           out = self.output_path
        else:
           out = os.path.dirname(os.path.abspath(self.output_path))

        #search for the chromosome length in the chr len files
        found = False
        try:
            #print '%s/../chrdesc/%s'%(os.path.dirname(__file__), self.species)
            for line in file('%s/../chrdesc/%s'%(os.path.dirname(__file__), self.species)):
                
                chrom, length = line.split()
                if chrom == chromosome:
                    found = True
                    self.chr_length = int(length)
        except IOError:
            pass #file not found, the warning will be printed
        
        if not found: self.logger.warning("The file containing %s length for assembly %s could not be found, an aproximation will be used"%(chromosome, self.species))
        self.result_log.write('---------------\n')
        self.result_log.write('%s\n'%(chromosome))
        self.result_log.write('---------------\n\n')
        self.result_log.write('Correction factor: %s\n\n'%(self.correction_factor))
        self.reads_per_bp =  self.total_bp_with_reads / self.chr_length*self.correction_factor
        
        p_nucleotide = 1.
        p_cluster = 1.
        p_numreads = 1.
        k = 0
        self.result_log.write('k\tcluster_length\tcluster_height\tnumreads\n')
        while ((self.absolute_max_numreads > k) or (self.absolute_max_height > k)) and k < self.height_limit:
            p_nucleotide -= Utils.poisson(k, self.reads_per_bp) #analisis nucleotide
            p_cluster -= Utils.poisson(k, self.acum_height/self.total_clusters) #analysis cluster
            p_numreads -= Utils.poisson(k, self.total_reads/self.total_clusters) #analysis numreads

            p_nucleotide = self._correct_bias(p_nucleotide)
            p_cluster = self._correct_bias(p_cluster)
            p_numreads = self._correct_bias(p_numreads)

            self.result_log.write('%s\t%.8f\t%.8f\t%.8f\n'%(k, p_nucleotide, p_cluster, p_numreads))
            
            if chromosome not in self.poisson_results['length'].keys() and p_nucleotide < self.p_value: #if we don't have a height k that is over the p_value yet, write it.
                self.poisson_results["length"][chromosome] = k

            if chromosome not in self.poisson_results['height'].keys() and p_cluster < self.p_value:
                self.poisson_results["height"][chromosome] = k
            
            if chromosome not in self.poisson_results['numreads'].keys() and p_numreads < self.p_value:
                self.poisson_results["numreads"][chromosome] = k

            if k not in self.maxheight_to_pvalue:
                self.maxheight_to_pvalue[k] = {}
            self.maxheight_to_pvalue[k][chromosome] = p_cluster

            if k not in self.numreads_to_pvalue:
                self.numreads_to_pvalue[k] = {}
            self.numreads_to_pvalue[k][chromosome] = p_numreads


            k+=1

    def poisson_retrieve_data(self, cluster):

        if self.estimate_frag_size: #We need to estimate the fragment size if not provided
            if cluster.tag_length:
                self.frag_size = cluster.tag_length
            else:
                self.frag_size = 100 
                
        acum_numreads = 0.
        self.total_clusters+=1
        for length, height in cluster:
            self.total_bp_with_reads+=length
            acum_numreads += length*height
        self.chr_length = max(self.chr_length, cluster.end)
        max_height = cluster.max_height()
        #numreads per cluster
        numreads_in_cluster = acum_numreads/self.frag_size
        self.total_reads += numreads_in_cluster
        self.absolute_max_numreads = max(numreads_in_cluster, self.absolute_max_numreads)
        #maxheight per cluster
        self.acum_height += max_height
        self.absolute_max_height = max(max_height, self.absolute_max_height)

    def process_cluster(self, cluster, output):
        if self.cluster._tag_cache:
            self.cluster._flush_tag_cache()

        if cluster.chromosome not in self.discarded_chromosomes and not cluster.is_empty():
            if self.previous_chr != cluster.chromosome: #A new chromosome has been detected
                if self.is_sorted and self.verbose:
                    print '%s...'%cluster.chromosome,
                    sys.stdout.flush()
                if self.do_poisson and not self.first_chr:
                    self.poisson_analysis(self.previous_chr)
                    self._init_poisson()
                self.previous_chr = cluster.chromosome
                if self.output_format == VARIABLE_WIG:
                    output.write('variableStep\tchrom=%s\tspan=%s\n'%(self.previous_chr, self.span))
                self.first_chr = False

            if self.do_subtract:
                cluster = self.subtract(cluster)

                for cluster in cluster.absolute_split(threshold=0):
                    self._post_subtract_process_cluster(cluster, output)
            else:
                self._post_subtract_process_cluster(cluster, output)



    def _post_subtract_process_cluster(self, cluster, output):
        if self.do_poisson:
            self.poisson_retrieve_data(cluster)

        if not (cluster.is_artifact() and DiscardArtifacts in self.operations) and not NoWrite in self.operations:
            if self.do_trim:
                cluster.trim(self.trim_absolute)

            if self.do_split:
                for subcluster in cluster.split(self.split_proportion, self.split_absolute): 
                    self.extract_and_write(subcluster, output)
            else:
                self.extract_and_write(cluster, output)


    def extract_and_write(self, cluster, output):
        """The line will be written to the file if the last conditions are met"""
        if not cluster.is_empty() and cluster.start > -1:
            output.write(cluster.write_line())

    def _current_directory(self):
        return os.path.abspath(os.path.dirname(self.current_output_path))

    def cut(self):
        """Discards the clusters that dont go past the threadshold calculated by the poisson analysis"""
        current_directory = self._current_directory()
        old_output = '%s/before_cut_%s'%(current_directory, os.path.basename(self.current_output_path))
        shutil.move(os.path.abspath(self.current_output_path), old_output)
        filtered_output = file(self.current_output_path, 'w+')
        if self.verbose:
            print "Filtering using",
            if self.poisson_test == 'height':
                print "cluster length..."
            else:
                print "number of reads per cluster..."

        unfiltered_output = file('%s/unfiltered_%s'%(current_directory, os.path.basename(self.current_output_path)), 'w+')
        if self.output_format == WIG or self.output_format == VARIABLE_WIG:
            wig_header = 'track type=wiggle_0\tname="%s"\tvisibility=full\n'%self.label
            filtered_output.write(wig_header)
            unfiltered_output.write(wig_header)
        cut_cluster = Cluster(read=self.output_format, write=self.output_format, rounding=self.rounding, read_half_open = self.open_output, write_half_open = self.open_output, tag_length=self.tag_length, span = self.span)
        print 'Writing filtered and unfiltered file...'
        for line in file(old_output):
            cut_cluster.clear()
            self.safe_read_line(cut_cluster, line)
            try:
                if self.poisson_test == 'height':
                    cut_cluster.p_value = self.maxheight_to_pvalue[int(round(cut_cluster.max_height()))][cut_cluster.chromosome]
                else:
                    cut_cluster.p_value = self.numreads_to_pvalue[int(round(cut_cluster.area()/self.frag_size))][cut_cluster.chromosome]
            except KeyError:
                cut_cluster.p_value = 0 #If the cluster is not in the dictionary, it means its too big, so the p_value will be 0

            try:
                if self.threshold:
                    thres = self.threshold
                else:
                    thres = self.poisson_results[self.poisson_test][cut_cluster.chromosome]
                    
                if cut_cluster.is_significant(thres):
                    filtered_output.write(cut_cluster.write_line())
                    
            except KeyError:
                pass

            if not cut_cluster.is_empty():
                unfiltered_output.write(cut_cluster.write_line()) 

        self._manage_temp_file(old_output)


    def __enrichment_A(self, region_a, total_reads_a, region_b, total_reads_b):
        return math.log(region_a.rpkm(total_reads_a)/region_b.rpkm(total_reads_b))


    def __enrichment_M(self, region_a, total_reads_a, region_b, total_reads_b):
        return (math.log(region_a.rpkm(total_reads_a))+math.log(region_b.rpkm(total_reads_b))/2)


    def enrichment(self):
        if self.verbose: print "Calculating enrichment in regions",
        file_a_reader = SortedFileClusterReader(self.current_experiment_path, self.experiment_format, cached=self.cached)
        file_b_reader = SortedFileClusterReader(self.current_control_path, self.experiment_format, cached=self.cached)
        if self.replica_a_path and self.replica_b_path:
            replica_a_reader = SortedFileClusterReader(self.current_replica_a_path, self.experiment_format, cached=self.cached)
            replica_b_reader = SortedFileClusterReader(self.current_replica_b_path, self.experiment_format, cached=self.cached)
            if self.verbose: print "using replicas..."
        else:
            if self.verbose: print "using swap..."
        self.total_reads_a = sum(1 for line in open(self.current_experiment_path))
        self.total_reads_b = sum(1 for line in open(self.current_control_path))     
        self.average_total_reads = (self.total_reads_a+self.total_reads_b)/2   
        real_A = []
        real_M = []
        swap_A = []
        swap_M = []
        for region_line in file(self.sorted_region_path):
            sregion = region_line.split()
            region_of_interest = Region(int(sregion[1]), int(sregion[2]), sregion[4], sregion[0])
            tags_a = file_a_reader.get_overlaping_clusters(region_of_interest, overlap=0.5)
            tags_b = file_b_reader.get_overlaping_clusters(region_of_interest, overlap=0.5)

            if tags_a or tags_b:
                for strand in (PLUS_STRAND, MINUS_STRAND):
                    region_a = region_of_interest.copy()
                    region_b = region_of_interest.copy()
                    region_a.add_tags(tags_a, strand=strand) #get only the tags of the desired strand
                    region_b.add_tags(tags_b, strand=strand) #get only the tags of the desired strand
                    swap1, swap2 = region_a.swap(region_b)
                    real_A.append(self.__enrichment_A(region_a, self.total_reads_a, region_b, self.total_reads_b))
                    real_M.append(self.__enrichment_M(region_a, self.total_reads_a, region_b, self.total_reads_b))
                    swap_A.append(self.__enrichment_A(swap1, self.average_total_reads, swap2, self.average_total_reads))
                    swap_M.append(self.__enrichment_M(swap1, self.average_total_reads, swap2, self.average_total_reads))
            
        try:
            from matplotlib.pyplot import plot, hist, show, legend
            plot(real_M, real_A, 'r.', label='Experiment')
            plot(swap_M, swap_A, 'b.', label='Swap')
            legend(loc=2)
            self._save_figure("enrichment_MA")
            hist(real_A, 100, color="r", histtype="step", label='Experiment')
            hist(swap_A, 100, color="b", histtype="step", label='Swap')
            legend(loc=2)
            self._save_figure("hist_D")
      
        except ImportError:
            print self.logger.warning('Pyicos can not find an installation of matplotlib, so no plot will be drawn for the strand correlation. If you want to get a plot with the correlation values, install the matplotlib library.')
            
    def _save_figure(self, figure_name):
        from matplotlib.pyplot import savefig, clf
        if os.path.dirname(self.current_output_path):
            figure_path = '%s/%s_%s.png'%(os.path.dirname(self.current_output_path), figure_name, os.path.basename(self.current_output_path))
        else:
            figure_path = '%s_%s.png'%(os.path.basename(figure_name), self.current_output_path)

        savefig(figure_path)
        clf()
        if self.verbose: print "%s figure saved to %s"%(figure_name, figure_path)


    def modfdr(self):
        print "\nRunning modfdr filter with %s p-value threshold and %s repeats..."%(self.p_value, self.repeats)
        old_output = '%s/before_modfdr_%s'%(self._current_directory(), os.path.basename(self.current_output_path))
        shutil.move(os.path.abspath(self.current_output_path), old_output)
        cluster_reader = SortedFileClusterReader(old_output, self.output_format, cached=self.cached)
        if self.masker_file:
            masker_reader = SortedFileClusterReader(self.masker_file, BED, cached=self.cached)
        filtered_output = file(self.current_output_path, 'w+')
        unfiltered_output = file('%s/unfiltered_%s'%(self._current_directory(), os.path.basename(self.current_output_path)), 'w+')
        for region_line in file(self.sorted_region_path):
            split_region_line = region_line.split()
            region = Region(split_region_line[1], split_region_line[2], chromosome=split_region_line[0])
            region.add_tags(cluster_reader.get_overlaping_clusters(region, overlap=0.000001), True)
            classified_clusters = region.get_FDR_clusters(self.repeats)
            for cluster in  classified_clusters:
                unfiltered_output.write(cluster.write_line())
                if cluster.p_value < self.p_value:
                    filtered_output.write(cluster.write_line())

        self._manage_temp_file(old_output)



    def strand_correlation(self):
        if self.verbose: print "Strand correlation analysis..."
        self.delta_results = dict()
        self.best_delta = -1
        positive_cluster = Cluster()
        negative_cluster = Cluster()
        positive_cluster_cache = [] #we are trying to hold to the previous cluster
        self.analyzed_pairs = 0.
        acum_length = 0.
        num_analyzed = 0
        for line in file(self.current_experiment_path):
            line_read = Cluster(read=self.experiment_format)
            line_read.read_line(line)
            if line_read.strand == '+':
                if  positive_cluster.intersects(line_read):
                     positive_cluster += line_read
                     acum_length += len(line_read)
                     num_analyzed += 1
                elif positive_cluster.is_empty() or positive_cluster.is_artifact():
                    positive_cluster = line_read.copy_cluster()
                elif not positive_cluster_cache:
                    positive_cluster_cache.append(line_read.copy_cluster())
                elif line_read.intersects(positive_cluster_cache[0]):
                    positive_cluster_cache[0] += line_read
                else:
                    positive_cluster_cache.append(line_read.copy_cluster())

            else:
                if negative_cluster.is_empty() or not negative_cluster.intersects(line_read):
                    if positive_cluster.max_height() > self.correlation_height and negative_cluster.max_height() > self.correlation_height: #if we have big clusters, correlate them
                        self._correlate_clusters(positive_cluster, negative_cluster)
                        if self.max_correlations < self.analyzed_pairs: #if we analyzed enough clusters, stop
                            break                     
                    negative_cluster = line_read.copy_cluster() #after correlating, select the next cluster
                else:
                    negative_cluster += line_read
                    acum_length += len(line_read)
                    num_analyzed += 1
                #advance in the positive cluster cache if its too far behind
                distance = negative_cluster.start-positive_cluster.start
                while distance > self.max_delta or positive_cluster.chromosome < negative_cluster.chromosome: # if the negative clusters are too far behind, empty the positive cluster
                    positive_cluster.clear()
                    if positive_cluster_cache:
                        positive_cluster = positive_cluster_cache.pop() 
                    else:
                        break

        #Use the results
        data = []
        max_delta = 0
        max_corr = -1
        average_len = acum_length/num_analyzed
        self.result_log.write("Strand Correlation\n")
        self.result_log.write("------------------\n\n")
        if self.verbose: print "Average analyzed length", average_len
        self.result_log.write("Average analyzed length:%s\n"%average_len)
        for delta in range(self.min_delta, self.max_delta, self.delta_step):
            if delta in self.delta_results:
                self.delta_results[delta]=self.delta_results[delta]/self.analyzed_pairs
                data.append(self.delta_results[delta])
                if self.delta_results[delta] > max_corr:
                    max_delta = delta
                    max_corr = self.delta_results[delta]
                #print 'Delta %s:%s'%(delta, self.delta_results[delta])
        if self.verbose: print 'Correlation test RESULT: You should extend this dataset to %s nucleotides'%(max_delta+average_len)
        self.result_log.write('Correlation test RESULT: You should extend this dataset to %s nucleotides\n'%(max_delta+average_len))
        self.frag_size = int(round(max_delta+average_len))
        if not data:
            if self.verbose: self.logger.warning('Not enough data to plot the correlation graph. Lower the threshold of the --height-filter flag')
        else: 
            try:
                import matplotlib.pyplot
                matplotlib.pyplot.plot(range(self.min_delta, self.max_delta), data)
                matplotlib.pyplot.plot()
                if os.path.dirname(self.current_output_path):
                    figure_path = '%s/%s.png'%(os.path.dirname(self.current_output_path), os.path.basename(self.current_output_path))
                else:
                    figure_path = '%s.png'%(os.path.basename(self.current_output_path))

                matplotlib.pyplot.savefig(figure_path)
                matplotlib.pyplot.clf()
                if self.verbose: print "Correlation figure saved to %s"%figure_path
                #matplotlib.pyplot.show()
                matplotlib.pyplot.clf()
            except ImportError:
                print self.logger.warning('Pyicos can not find an installation of matplotlib, so no plot will be drawn for the strand correlation. If you want to get a plot with the correlation values, install the matplotlib library.')


    def _correlate_clusters(self, positive_cluster, negative_cluster):
        distance = negative_cluster.end-positive_cluster.start
        if (distance < self.max_delta and distance > self.min_delta) and (positive_cluster.chromosome == negative_cluster.chromosome):
            self.analyzed_pairs+=1
            for delta in range(self.min_delta, self.max_delta+1, self.delta_step):
                r_squared = self._analize_paired_clusters(positive_cluster, negative_cluster, delta)**2
                if delta not in self.delta_results:
                    self.delta_results[delta] = r_squared
                else:
                    self.delta_results[delta] += r_squared


    def _analize_paired_clusters(self, positive_cluster, negative_cluster, delta):
        #from scipy.stats.stats import pearsonr Abandoned scipy
        positive_array = []
        negative_array = [] 
        #delta correction
        corrected_positive_start = positive_cluster.start + delta
        #add zeros at the start of the earliest cluster
        if corrected_positive_start > negative_cluster.start:
            self.__add_zeros(positive_array, corrected_positive_start - negative_cluster.start)
        elif negative_cluster.start > corrected_positive_start:  
            self.__add_zeros(negative_array, negative_cluster.start - corrected_positive_start)
        #add the values of the clusters
        positive_array.extend(positive_cluster.get_heights())
        negative_array.extend(negative_cluster.get_heights())
        #add the zeros at the end of the shortest array
        if len(positive_array) > len(negative_array):
            self.__add_zeros(negative_array, len(positive_array) - len(negative_array))
        elif len(positive_array) < len(negative_array):
            self.__add_zeros(positive_array, len(negative_array) - len(positive_array))
        return Utils.pearson(negative_array, positive_array)

    
    def __add_zeros(self, array, num_zeros):
        for i in range(0, num_zeros):
            array.append(0)



