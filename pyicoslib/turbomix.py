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
from datetime import datetime
from tempfile import gettempdir

#pyicos stuff
from core import Cluster, Region, InvalidLine, InsufficientData, ConversionNotSupported
from defaults import *
import utils

class OperationFailed(Exception):
    pass                 

class Turbomix:
    """
    This class is the pipeline that makes possible the different combination of operations. 
    It has different switches that are activated by the list 'self.operations'.
    """
    logging_format= "%(asctime)s - %(levelname)s - %(message)s"
    logging.basicConfig(filename="pyicos.log", format=logging_format)
    logger = logging.getLogger("pyicos.log")
    logger.setLevel(logging.WARNING)
    ch = logging.StreamHandler()
    ch.setLevel(logging.WARNING)
    formatter = logging.Formatter("%(levelname)s - %(message)s")
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    invalid_count = 0
    invalid_limit = 2000
    
    def __init__(self, experiment_path, output_path, experiment_format=BED, output_format=PK, label=LABEL, 
                 open_experiment=OPEN_EXPERIMENT, open_output=OPEN_OUTPUT, debug = DEBUG, rounding=ROUNDING, tag_length = TAG_LENGTH, discarded_chromosomes = REMLABELS,
                 control_path = CONTROL, control_format = PK, open_control = OPEN_CONTROL, region_path = REGION, region_format = PK, 
                 open_region=OPEN_REGION, span = SPAN, frag_size = FRAG_SIZE, p_value = P_VALUE, height_limit = HEIGHT_LIMIT, correction_factor = CORRECTION,
                 trim_percentage=TRIM_PROPORTION, no_sort=NO_SORT, duplicates=DUPLICATES, threshold=THRESHOLD, trim_absolute=TRIM_ABSOLUTE, max_delta=MAX_DELTA, 
                 min_delta=MIN_DELTA, correlation_height=HEIGHT_FILTER, delta_step=DELTA_STEP, verbose=VERBOSE, species=SPECIES, cached=CACHED, split_proportion=SPLIT_PROPORTION, 
                 split_absolute=SPLIT_ABSOLUTE, repeats=REPEATS, masker_file=MASKER_FILE, max_correlations=MAX_CORRELATIONS, keep_temp=KEEP_TEMP, experiment_b_path=EXPERIMENT,
                 replica_a_path=EXPERIMENT, replica_b_path=EXPERIMENT, poisson_test=POISSONTEST, stranded_analysis=STRANDED_ANALYSIS, proximity=PROXIMITY, 
                 postscript=POSTSCRIPT, showplots=SHOWPLOTS, plot_path=PLOT_PATH, no_pseudocount=NOPSEUDOCOUNT, simple_counts=SIMPLECOUNTS, label1=LABEL1, 
                 label2=LABEL2, numbins=NUMBINS, zscore=ZSCORE):

        self.__dict__.update(locals())

        self.is_sorted = False
        self.temp_experiment = False #Indicates if temp files where created for the experiment
        self.temp_control = False #Indicates if temporary files where created for the control
        self.temp_replica_a = False 
        self.temp_replica_b = False 
        self.operations = []
        self.previous_chr = None
        self.open_region = False
        self.safe_reader = utils.SafeReader(self.verbose)        

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
            utils.list_available_formats()
        #duplicates flag
        self.previous_start = 0
        self.previous_end = 0
        self.previous_chr = ''
        self.duplicates_found = 0
        #poisson stuff
        self.first_chr = True
        self._init_poisson()
        self.poisson_results = {'length': defaultdict(), 'height': defaultdict(), 'numtags': defaultdict()}
        self.maxheight_to_pvalue = {}
        self.numtags_to_pvalue = {}
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
        if FILTER in self.operations and POISSON not in self.operations and not self.threshold:
            print "Can't do Filter without Poisson or a fixed threshold\n"
            sys.exit(0)
        elif (SUBTRACT or SPLIT or POISSON) in self.operations and (self.output_format not in CLUSTER_FORMATS):
            print "Can't get the output as tag format (eland, bed) for Subtract, Split, Poisson filtering please use a clustered format %s\n"%CLUSTER_FORMATS
            sys.exit(0)
        elif EXTEND in self.operations and self.experiment_format in CLUSTER_FORMATS:
            print "Can't extend if the experiment is a clustered format ",
            print CLUSTER_FORMATS
            sys.exit(0)
        elif STRAND_CORRELATION in self.operations and self.experiment_format in CLUSTER_FORMATS:
            print "Can't perform strand correlation operation if the experiment is a clustered format ",
            print CLUSTER_FORMATS
            sys.exit(0)
        
        elif (ENRICHMENT and ModFDR) in self.operations:
            print
            print "Enrichment and ModFDR operations are not compatible"
            sys.exit(0)

        elif (ENRICHMENT and POISSON) in self.operations:
            print
            print "Enrichment and Poisson operations are not compatible"
            sys.exit(0)             
    
        
    def _add_slash_to_path(self, path):
        return utils.add_slash_to_path(path)

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
        
            if self.experiment_format != self.output_format and not NOWRITE in self.operations:
                print 'Convert from %s to %s.'%(self.experiment_format, self.output_format)

        if STRAND_CORRELATION in self.operations:
            self.result_log.write('Strand correlation: %.0f pairs analyzed, estimated a %s nucleotides extension value'%(self.analyzed_pairs, self.frag_size))
            if self.verbose: print 'Strand correlation'

        if EXTEND in self.operations:
            self.result_log.write('EXTEND to %s\n'%self.frag_size)
            if self.verbose: print 'EXTEND to %s'%self.frag_size

        if self.do_subtract:
            self.result_log.write('SUBTRACT\n')
            if self.verbose: print 'SUBTRACT'

        if self.discarded_chromosomes:
            self.result_log.write('Discard tags: %s\n'%self.discarded_chromosomes)
            if self.verbose: print 'Discard tags: %s'%self.discarded_chromosomes

        if self.do_heuremove:
            self.result_log.write('Heuristic Remove from %s\n'%self.region_path)
            if self.verbose: print 'Heuristic Remove from %s'%self.region_path

        if self.do_split:
            self.result_log.write('SPLIT\n')
            if self.verbose: print 'SPLIT'

        if DISCARD_ARTIFACTS in self.operations:
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
        if not NOWRITE in self.operations and self.verbose:
            print 'Output at: %s'%(output_path)

    def start_operation_message(self):
        if self.verbose:
            if self.experiment_format != self.output_format and not NOWRITE in self.operations:
                print 'Converting file %s from %s to %s...'%(self.current_experiment_path, self.experiment_format, self.output_format)
            else:
                print 'Reading file %s as a %s file...'%(self.current_experiment_path, self.experiment_format)

            if self.current_control_path:
                print 'Control file:%s'%self.current_control_path
                if self.do_normalize:
                    print 'The file %s will be normalized to match %s'%(self.current_experiment_path, self.current_control_path)
                if self.do_subtract:
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
        self.safe_reader.safe_read_line(cluster, line)

    def run(self):
        self.do_subtract = (SUBTRACT in self.operations and self.control_path)
        self.do_normalize = (NORMALIZE in self.operations and self.control_path)
        self.do_heuremove = (REMOVE_REGION in self.operations and self.region_path)
        self.do_poisson = POISSON in self.operations
        self.do_split = SPLIT in self.operations
        self.do_trim = TRIM in self.operations
        self.do_cut = FILTER in self.operations
        self.do_extend = EXTEND in self.operations and not self.sorted_by_pyicos #If the experiment was sorted by Pyicos, it was already extended before, so dont do it again
        self.do_discard = DISCARD_ARTIFACTS in self.operations
        self.do_dupremove = REMOVE_DUPLICATES in self.operations
        self.tag_to_cluster = (not self.experiment_format in CLUSTER_FORMATS and self.output_format in CLUSTER_FORMATS)
        if self.verbose:
            print '\n\nPyicos running...'
        if self.control_path:
            self.process_all_files_paired(self.experiment_path, self.control_path)
        elif self.experiment_b_path:
            self.process_all_files_paired(self.experiment_path, self.experiment_b_path)
        elif self.plot_path:
            self.process_all_files_recursive(self.plot_path, self.output_path)
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
        #TODO normalize should go here too
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
        if self.tag_to_cluster or self.do_subtract or self.do_heuremove or self.do_dupremove or ModFDR in self.operations or ENRICHMENT in self.operations or REMOVE_REGION in self.operations:
            self.experiment_preprocessor = utils.BigSort(self.experiment_format, self.open_experiment, self.frag_size, 'experiment', verbose=self.verbose)
            self.experiment_b_preprocessor = utils.BigSort(self.experiment_format, self.open_experiment, self.frag_size, 'experiment_b', verbose=self.verbose)
            self.replica_a_preprocessor = utils.BigSort(self.experiment_format, self.open_experiment, self.frag_size, 'replica_a', verbose=self.verbose)
            self.replica_b_preprocessor = utils.BigSort(self.experiment_format, self.open_control, self.frag_size, 'replica_b', verbose=self.verbose)
            self.control_preprocessor = utils.BigSort(self.control_format, self.open_control, self.frag_size, 'control', verbose=self.verbose)

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
            
            if ENRICHMENT in self.operations:  
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
            self.region_preprocessor = utils.BigSort(self.region_format, self.open_region, None, 'region', verbose=self.verbose)
            self.sorted_region_file = self.region_preprocessor.sort(self.region_path, None, self.get_lambda_func(BED))
            self.sorted_region_path = self.sorted_region_file.name

    def add_operation(self, operation):
        create_operation(operation)

    def operate(self, experiment_path, control_path=None, output_path=None):
        """Operate expects single paths, not directories. Its called from run() several times if the experiment for picos is a directory"""
        try:
            self.i_cant_do()
            #per operation variables
            self.previous_chr = None
            self.real_experiment_path = experiment_path #for refering to the files after all the intermediate files
            self.real_control_path = control_path #for refering to the files after all the intermediate files
            self.current_experiment_path = experiment_path
            self.current_control_path = control_path
            self.current_output_path = output_path
            self.cluster.clear()
            self.cluster_aux.clear()
            self.result_log = file('%s/pyicos_report_%s.txt'%(self._output_dir(), os.path.basename(self.current_output_path)), 'wb')
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
                self.control_reader = utils.SortedFileClusterReader(self.current_control_path, self.control_format, rounding=self.rounding, cached=self.cached)
            
            if self.region_path:
                self.region_reader = utils.SortedFileClusterReader(self.sorted_region_path, BED, rounding=self.rounding, cached=self.cached)

            if STRAND_CORRELATION in self.operations:
                self.strand_correlation()
            
            if self.do_dupremove:
                self.filter_files()

            if NORMALIZE in self.operations:
                self.normalize()

            if self.do_cut: #if we cut, we will round later
                self.cluster.rounding = False
                self.cluster_aux.rounding = False

            if ((not NOWRITE in self.operations) or (NOWRITE in self.operations and POISSON in self.operations)) and not ENRICHMENT in self.operations: 
                self.process_file()
            
            if self.do_poisson: #extract info from the last chromosome and print the thresholds
                self.poisson_analysis(self.previous_chr)
                print '\nCluster threadsholds for p-value %s:'%self.p_value
                self.result_log.write('\nCluster threadsholds for p-value %s:\n'%self.p_value)
                for chromosome, k in self.poisson_results[self.poisson_test].items():
                    print '%s: %s'%(chromosome, k)
                    self.result_log.write('%s: %s\n'%(chromosome, k))


            #Mutually exclusive final operations
            if ENRICHMENT in self.operations:
                self.plot_path = self.enrichment()

            elif self.do_cut: 
                self.cut()

            elif ModFDR in self.operations:
                self.modfdr()

            #Plot and summary operations
            if PLOT in self.operations:
                  if self.verbose: print "Plotting..."
                  self.plot_enrichment(self.plot_path)              

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
        #load the first read
        while self.cluster.is_empty():
            self.cluster_aux2.clear()
            self.read_and_preprocess(self.cluster_aux2, experiment.next())
            if not self.cluster_aux2.is_empty():
                self.cluster = self.cluster_aux2.copy_cluster()
            
        for line in experiment:
            self.cluster_aux.clear()
            self.read_and_preprocess(self.cluster_aux, line)
            if not self.cluster_aux.is_empty():
                if self.cluster.intersects(self.cluster_aux) or self.cluster_aux.is_contiguous(self.cluster) or self.cluster.is_contiguous(self.cluster_aux):
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
        if NOWRITE in self.operations:
            output = None
        else:
            output = file(self.current_output_path, 'wb')
        experiment = file(self.current_experiment_path, 'rb')

        if self.output_format == WIG or self.output_format == VARIABLE_WIG:
            output.write('track type=wiggle_0\tname="%s"\tvisibility=full\n'%self.label)

        if self.output_format in CLUSTER_FORMATS and not ModFDR in self.operations and not ENRICHMENT in self.operations:
            if not self.experiment_format in CLUSTER_FORMATS and self.verbose:
                print 'Clustering reads...'
            self._to_cluster_conversion(experiment, output)
        else:
            self._to_read_conversion(experiment, output)

        if not NOWRITE in self.operations:
            output.flush()
            output.close()

    def subtract(self, cluster):
        region = Region(cluster.start, cluster.end, cluster.chromosome)
        over = self.control_reader.get_overlaping_clusters(region, overlap=0.0000001)
        for tag in over:
            if self.do_extend:
                tag.extend(self.frag_size)

       
        region.add_tags(over, True)
        meta = region.get_metacluster()
        cluster -= meta
        return cluster

    def remove_regions(self, cluster):
        region = Region(cluster.start, cluster.end, cluster.chromosome)
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
        self.absolute_max_numtags = 0
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
        p_numtags = 1.
        k = 0
        self.result_log.write('k\tcluster_length\tcluster_height\tnumtags\n')
        while ((self.absolute_max_numtags > k) or (self.absolute_max_height > k)) and k < self.height_limit:
            p_nucleotide -= utils.poisson(k, self.reads_per_bp) #analisis nucleotide
            p_cluster -= utils.poisson(k, self.acum_height/self.total_clusters) #analysis cluster
            p_numtags -= utils.poisson(k, self.total_reads/self.total_clusters) #analysis numtags

            p_nucleotide = self._correct_bias(p_nucleotide)
            p_cluster = self._correct_bias(p_cluster)
            p_numtags = self._correct_bias(p_numtags)

            self.result_log.write('%s\t%.8f\t%.8f\t%.8f\n'%(k, p_nucleotide, p_cluster, p_numtags))
            
            if chromosome not in self.poisson_results['length'].keys() and p_nucleotide < self.p_value: #if we don't have a height k that is over the p_value yet, write it.
                self.poisson_results["length"][chromosome] = k

            if chromosome not in self.poisson_results['height'].keys() and p_cluster < self.p_value:
                self.poisson_results["height"][chromosome] = k
            
            if chromosome not in self.poisson_results['numtags'].keys() and p_numtags < self.p_value:
                self.poisson_results["numtags"][chromosome] = k

            if k not in self.maxheight_to_pvalue:
                self.maxheight_to_pvalue[k] = {}
            self.maxheight_to_pvalue[k][chromosome] = p_cluster

            if k not in self.numtags_to_pvalue:
                self.numtags_to_pvalue[k] = {}
            self.numtags_to_pvalue[k][chromosome] = p_numtags


            k+=1

    def poisson_retrieve_data(self, cluster):

        if self.estimate_frag_size: #We need to estimate the fragment size if not provided
            if cluster.tag_length:
                self.frag_size = cluster.tag_length
            else:
                self.frag_size = 100 
                
        acum_numtags = 0.
        self.total_clusters+=1
        for length, height in cluster:
            self.total_bp_with_reads+=length
            acum_numtags += length*height
        self.chr_length = max(self.chr_length, cluster.end)
        max_height = cluster.max_height()
        #numtags per cluster
        numtags_in_cluster = acum_numtags/self.frag_size
        self.total_reads += numtags_in_cluster
        self.absolute_max_numtags = max(numtags_in_cluster, self.absolute_max_numtags)
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

        if not (cluster.is_artifact() and DISCARD_ARTIFACTS in self.operations) and not NOWRITE in self.operations:
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
                print "cluster height..."
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
                    cut_cluster.p_value = self.numtags_to_pvalue[int(round(cut_cluster.area()/self.frag_size))][cut_cluster.chromosome]
            except KeyError:
                cut_cluster.p_value = 0 #If the cluster is not in the dictionary, it means its too big, so the p_value will be 0

            try:
                if self.threshold:
                    thres = self.threshold
                else:
                    thres = self.poisson_results[self.poisson_test][cut_cluster.chromosome]
                    
                if cut_cluster.is_significant(thres, self.poisson_test):
                    filtered_output.write(cut_cluster.write_line())
                    
            except KeyError:
                pass

            if not cut_cluster.is_empty():
                unfiltered_output.write(cut_cluster.write_line()) 

        self._manage_temp_file(old_output)

    def _output_dir(self):
        """Returns the output directory"""
        path = os.path.dirname(os.path.realpath(self.current_output_path))
        if not os.path.exists(path):
            os.makedirs(path)
        return path

    def _region_from_sline(self, sline):
        if self.stranded_analysis:
            return Region(int(sline[1]), int(sline[2]), chromosome=sline[0], strand=sline[5])
        else:
            return Region(int(sline[1]), int(sline[2]), chromosome=sline[0]) 
     

    def calculate_region(self):
        """
        Calculate a region file using the reads present in the both main files to analyze. 
        """
        #TODO stranded support
        if self.verbose: print 'Generating region file...'
        self.sorted_region_path = '%s/calcregion_%s.txt'%(self._output_dir(), os.path.basename(self.current_output_path))
        region_file = open(self.sorted_region_path, 'wb')
        dual_reader = utils.DualSortedReader(self.current_experiment_path, self.current_control_path, self.experiment_format, self.verbose) 

        self.calculate_region_notstranded(dual_reader, region_file)    

    def calculate_region_notstranded(self, dual_reader, region_file):
        calculated_region = Region()        
        for line in dual_reader:
            if not calculated_region: #first region only
                calculated_region = self._region_from_sline(line.split())
                calculated_region.end += self.proximity
            else:
                new_region = self._region_from_sline(line.split())
                new_region.end += self.proximity
                if calculated_region.overlap(new_region):
                    calculated_region.join(new_region)
                else: 
                    calculated_region.end -= self.proximity
                    region_file.write(calculated_region.write())                         
                    calculated_region = new_region.copy()

    def calculate_region_stranded(self, dual_reader, region_file):
        region_file = open(self.sorted_region_path, 'wb')
        calculated_region = Region()        
        dual_reader = utils.DualSortedReader(self.current_experiment_path, self.current_control_path, self.experiment_format, self.verbose)     
        for line in dual_reader:
            if not calculated_region: #first region only
                calculated_region = self._region_from_sline(line.split())
                calculated_region.end += self.proximity
            else:
                new_region = self._region_from_sline(line.split())
                new_region.end += self.proximity
                if calculated_region.overlap(new_region):
                    calculated_region.join(new_region)
                else: 
                    calculated_region.end -= self.proximity
                    region_file.write(calculated_region.write())                         
                    calculated_region = new_region.copy()           


    def plot_enrichment(self, file_path):
        try:
            if self.postscript:
                import matplotlib
                matplotlib.use("PS")

            from matplotlib.pyplot import plot, hist, show, legend, figure, xlabel, ylabel
            from matplotlib import rcParams
            rcParams['legend.fontsize'] = 8
            #decide labels
            if self.label1:
                label_main = self.label1
            else:
                if self.real_control_path and self.real_experiment_path:
                    label_main = '%s VS %s'%(os.path.basename(self.real_experiment_path), os.path.basename(self.real_control_path))
                else:
                    label_main = "A VS B"

            if self.label2:
                label_control = self.label2
            else:                
                if self.replica_a_path:
                    label_control = '%s(A) VS %s(A)'%(os.path.basename(self.real_experiment_path), os.path.basename(self.real_control_path))
                else:
                    label_control = 'Swap1 VS Swap2'

            A = []
            A_prime = []
            M = []
            M_significant = []
            A_significant = []
            M_prime = []
            figure(figsize=(8,6))

            for line in open(file_path): 
                sline = line.split()
                if abs(float(sline[16])) < self.zscore:
                    M.append(float(sline[5]))                
                    A.append(float(sline[4]))
                else:
                    M_significant.append(float(sline[5]))                
                    A_significant.append(float(sline[4]))     
               
                M_prime.append(float(sline[11]))    
                A_prime.append(float(sline[10]))
 
            plot(A, M, 'y.', label=label_main)
            plot(A_significant, M_significant, 'r.', label="%s (significant)"%label_main)
            plot(A_prime, M_prime, 'b.', label=label_control)
            xlabel('A')
            ylabel('M')
            legend(bbox_to_anchor=(0., 1.01, 1., .101), loc=3, ncol=1, mode="expand", borderaxespad=0.)
            self._save_figure("enrichment_MA")

            M.extend(M_significant)
            hist(M, 100, color="r", histtype="step", label=label_main)
            hist(M_prime, 100, color="b", histtype="step", label=label_control)
            xlabel('M')
            ylabel('')
            legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3, ncol=1, mode="expand", borderaxespad=0.)
            self._save_figure("hist_A")

        except ImportError:
            print self.logger.warning('Pyicos can not find an installation of matplotlib, so no plot will be drawn. If you want to get a plot with the correlation values, install the matplotlib library.')

    def __enrichment_M(self, rpkm_a, rpkm_b):
        return math.log(rpkm_a/float(rpkm_b), 2)

    def __enrichment_A(self, rpkm_a, rpkm_b):
        return (math.log(float(rpkm_a), 2)+math.log(float(rpkm_b), 2))/2

    def enrichment(self):
        use_pseudocount = not self.no_pseudocount
        use_replica = bool(self.replica_a_path)

        if self.verbose: print "Calculating enrichment in regions",
        file_a_reader = utils.SortedFileClusterReader(self.current_experiment_path, self.experiment_format, cached=self.cached)
        file_b_reader = utils.SortedFileClusterReader(self.current_control_path, self.experiment_format, cached=self.cached)
        out_file = open(self.current_output_path, 'wb')
        if use_replica:
            replica_a_reader = utils.SortedFileClusterReader(self.current_replica_a_path, self.experiment_format, cached=self.cached)
            if self.verbose: print "using replicas..."
        else:
            if self.verbose: print "using swap..."

        if self.sorted_region_path:
            print 'Using region file %s'%self.region_path
        else:
            self.calculate_region() #create region file semi automatically

        if self.verbose: print "... counting number of lines in files..."
        self.total_reads_a = sum(1 for line in open(self.current_experiment_path))
        self.total_reads_b = sum(1 for line in open(self.current_control_path))
        if use_replica:
            self.total_reads_replica_a = sum(1 for line in open(self.replica_a_path))
        self.total_regions = sum(1 for line in open(self.sorted_region_path))
        self.average_total_reads = (self.total_reads_a+self.total_reads_b)/2
        enrichment_result = [] #This will hold the chromosome, start and end of the region, plus the A, M, 'M and 'A
        
        regions_analyzed_count = 0
        if self.verbose: print "... analyzing regions..."
        for region_line in file(self.sorted_region_path):
            region_of_interest = self._region_from_sline(region_line.split())
            tags_a = file_a_reader.get_overlaping_clusters(region_of_interest, overlap=0.5)
            tags_b = file_b_reader.get_overlaping_clusters(region_of_interest, overlap=0.5)
            replica_a = None
            replica_tags = None
            swap1 = None
            swap2 = None
            if use_replica:            
                replica_tags = replica_a_reader.get_overlaping_clusters(region_of_interest, overlap=0.5)

            #if we are using pseudocounts, use the union, use the intersection otherwise
            if (use_pseudocount and (tags_a or tags_b)) or (not use_pseudocount and tags_a and tags_b): 
                region_a = region_of_interest.copy()
                region_b = region_of_interest.copy()
                region_a.add_tags(tags_a)
                region_b.add_tags(tags_b)
                if self.simple_counts:
                    rpkm_a = region_a.numtags(use_pseudocount)/self.total_reads_a
                    rpkm_b = region_b.numtags(use_pseudocount)/self.total_reads_b
                else:
                    rpkm_a = region_a.rpkm(self.total_reads_a, self.total_regions, use_pseudocount)
                    rpkm_b = region_b.rpkm(self.total_reads_b, self.total_regions, use_pseudocount)    


                A = self.__enrichment_A(rpkm_a, rpkm_b)
                M = self.__enrichment_M(rpkm_a, rpkm_b)  
                if use_replica:
                    replica_a = region_of_interest.copy()
                    replica_a.add_tags(replica_tags)
                    count_1 = len(tags_a)
                    count_2 = len(replica_tags)
                    total_1 = self.total_reads_a
                    total_2 = self.total_reads_replica_a
                    region_rpkm = rpkm_a
                    if self.simple_counts:
                        replica_a.numtags(use_pseudocount)
                    else:
                        replica_rpkm = replica_a.rpkm(self.total_reads_replica_a, self.total_regions, use_pseudocount)
                else:
                    swap1, swap2 = region_a.swap(region_b)
                    count_1 = len(swap1.tags)
                    count_2 = len(swap2.tags)
                    total_1 = total_2 = self.average_total_reads
                    if self.simple_counts: 
                        region_rpkm = swap1.numtags(use_pseudocount)/self.average_total_reads
                        replica_rpkm = swap2.numtags(use_pseudocount)/self.average_total_reads                     
                    else:                  
                        region_rpkm = swap1.rpkm(self.average_total_reads, self.total_regions, use_pseudocount)
                        replica_rpkm = swap2.rpkm(self.average_total_reads, self.total_regions, use_pseudocount)

                #if there is no data in the replica or in the swap and we are not using pseudocounts, dont write the data 
                if use_pseudocount or (use_replica and replica_tags) or (not use_replica and swap1.tags and swap2.tags):
                    M_prime = self.__enrichment_M(region_rpkm, replica_rpkm)
                    A_prime = self.__enrichment_A(region_rpkm, replica_rpkm)             
                    enrichment_result.append({'chr':region_of_interest.chromosome, 'start':region_of_interest.start, 'end':region_of_interest.end, 'strand': region_of_interest.strand,
                                              'A':A, 'M':M, 'total_A': self.total_reads_a, 'total_b':self.total_reads_b, 'num_tags_a':len(tags_a), 
                                              'num_tags_b':len(tags_b), 'A_prime':A_prime, 'M_prime':M_prime, 'total_1':total_1, 'total_2':total_2, 
                                              'count_1':count_1, 'count_2':count_2, 'z_score': 0})
                    regions_analyzed_count += 1


        points = self.calculate_zscore(enrichment_result)
        self.__sub_enrich_write(enrichment_result, out_file)
        out_file.flush()
        out_file.close()
        if self.verbose: 
            print "%s regions analyzed.\nEnrichment result saved to %s"%(regions_analyzed_count, self.current_output_path)
        return out_file.name

    def calculate_zscore(self, enrichment_result):
        if self.verbose: "... calculating z-score..."
        bin_size = int(len(enrichment_result)/(self.numbins+1)) #Make the bins proportionals #+1, since we are going to fuse the last 2
        enrichment_result.sort(key=lambda x:(x["A_prime"]))
        points = []
        #get the standard deviations
        for i in range(0, len(enrichment_result)-bin_size, bin_size):
            #get the slice
            if i+2*bin_size < len(enrichment_result):
                result_chunk = enrichment_result[i:i+bin_size]  
            else:
                result_chunk = enrichment_result[i:] #fuse the last two chunks together so there is not a small minichunk 

            #retrieve the values
            mean_acum = 0
            Ms_replica = []
            for entry in result_chunk:
                mean_acum += entry["M_prime"]
                Ms_replica.append(entry["M_prime"])

            #add them to the points of mean and sd
            mean = mean_acum/len(result_chunk)
            #sd = (sum((x - mean)**2 for x in Ms_replica))/len(result_chunk) #too fancy 
            sd_acum = 0  
            for x in Ms_replica:
                sd_acum += (x - mean)**2

            sd = sd_acum/len(result_chunk)
            points.append([result_chunk[-1]["A_prime"], mean, sd]) #The maximum A of the chunk, the mean and the standard deviation     
            print result_chunk[-1]["A_prime"], mean, sd

        #update z scores
        for entry in enrichment_result:
            for i in range(0, len(points)):
                if entry["A"] < points[i][0]:
                    self.__sub_zscore(entry, points, i)
                    continue #found it, leave go to the next
               
                self.__sub_zscore(entry, points, -1) #the value is the biggest, use the last break

        enrichment_result.sort(key=lambda x:(x["chr"], x["start"], x["end"]))
        #return points
        
    def __sub_zscore(self, entry, points, i):
        if points[i][2] > 0:
            entry["z_score"] = ((entry["M"]-points[i][1])/points[i][2]) 
        else:
            entry["z_score"] = 0 #This points are weird anyway

    def __sub_enrich_write(self, er, out_file):
        for d in er:
            out_file.write("%s\n"%"\t".join([d['chr'], str(d['start']), str(d['end']), str(d['strand']), str(d['A']), str(d['M']), str(d['total_A']), str(d['total_b']), str(d['num_tags_a']), 
                                            str(d['num_tags_b']), str(d['A_prime']), str(d['M_prime']), str(d['total_1']), str(d['total_2']), str(d['count_1']), 
                                            str(d['count_2']), str(d['z_score'])]))

    def _save_figure(self, figure_name):
        if self.postscript:
            exten = 'ps'
        else:
            exten = 'png'
        from matplotlib.pyplot import savefig, clf, show
        if self.showplots:
            show()
        else:
            figure_path = '%s/%s_%s.%s'%(self._output_dir(), figure_name, os.path.basename(self.current_output_path), exten)
            savefig(figure_path)
            if self.verbose: print "%s figure saved to %s"%(figure_name, figure_path)
        clf()
        


    def modfdr(self):
        print "\nRunning modfdr filter with %s p-value threshold and %s repeats..."%(self.p_value, self.repeats)
        old_output = '%s/before_modfdr_%s'%(self._current_directory(), os.path.basename(self.current_output_path))
        shutil.move(os.path.abspath(self.current_output_path), old_output)
        cluster_reader = utils.SortedFileClusterReader(old_output, self.output_format, cached=self.cached)
        if self.masker_file:
            masker_reader = utils.SortedFileClusterReader(self.masker_file, BED, cached=self.cached)
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
        positive_cluster = Cluster(cached=True)
        negative_cluster = Cluster(cached=True)
        positive_cluster_cache = [] #we are trying to hold to the previous cluster
        self.analyzed_pairs = 0.
        acum_length = 0.
        num_analyzed = 0
        for line in file(self.current_experiment_path):
            line_read = Cluster(read=self.experiment_format)
            line_read.read_line(line)
            if line_read.strand == PLUS_STRAND:
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
        if num_analyzed:
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
                if self.postscript:
                    import matplotlib
                    matplotlib.use("PS")

                import matplotlib.pyplot
                matplotlib.pyplot.plot(range(self.min_delta, self.max_delta), data)
                matplotlib.pyplot.plot()
                plt.xlabel("Shift length between + and - clusters")
                plt.ylabel("Pearson correlation coefficient")
                if os.path.dirname(self.current_output_path):
                    figure_path = '%s/%s.%s'%(os.path.dirname(self.current_output_path), os.path.basename(self.current_output_path), ext)
                else:
                    figure_path = '%s.%s'%(os.path.basename(self.current_output_path), ext)

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
        return utils.pearson(negative_array, positive_array)

    
    def __add_zeros(self, array, num_zeros):
        for i in range(0, num_zeros):
            array.append(0)



