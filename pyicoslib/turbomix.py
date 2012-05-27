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
import bam

try:
    from shutil import move
except:
    print "WARNING: Couldn't import shutil, using os.rename instead"
    from os import rename as move 



class OperationFailed(Exception):
    pass                 

class Turbomix:
    """
    This class is the pipeline that makes possible the different combination of operations. 
    It has different switches that are activated by the list 'self.operations'.
    """
    
    def set_logger(self, verbose, debug):
        logging_format= "%(asctime)s (PID:%(process)s) - %(levelname)s - %(message)s"
        logging.basicConfig(filename="pyicos.log", format=logging_format)
        self.logger = logging.getLogger("pyicos.log")
        if self.debug:
            self.logger.setLevel(logging.DEBUG)
        else:
            self.logger.setLevel(logging.INFO)

        ch = logging.StreamHandler()
        if self.debug: 
            ch.setLevel(logging.DEBUG)
        elif self.verbose:
            ch.setLevel(logging.INFO)
        else:
            ch.setLevel(logging.WARNING)

        formatter = logging.Formatter("%(asctime)s (PID:%(process)s) - %(levelname)s - %(message)s")
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)


    def __init__(self, experiment_path, output_path, experiment_format=BED, output_format=PK, label=LABEL, 
                 open_experiment=OPEN_EXPERIMENT, open_output=OPEN_OUTPUT, debug = DEBUG, rounding=ROUNDING, tag_length = TAG_LENGTH, discarded_names = REMLABELS,
                 control_path = CONTROL, control_format = PK, open_control = OPEN_CONTROL, region_path = REGION, region_format = PK, 
                 open_region=OPEN_REGION, span = SPAN, frag_size = FRAG_SIZE, p_value = P_VALUE, height_limit = HEIGHT_LIMIT, correction_factor = CORRECTION,
                 trim_percentage=TRIM_PROPORTION, no_sort=NO_SORT, duplicates=DUPLICATES, threshold=THRESHOLD, trim_absolute=TRIM_ABSOLUTE, max_delta=MAX_DELTA, 
                 min_delta=MIN_DELTA, correlation_height=HEIGHT_FILTER, delta_step=DELTA_STEP, verbose=VERBOSE, species=SPECIES, cached=CACHED, split_proportion=SPLIT_PROPORTION, 
                 split_absolute=SPLIT_ABSOLUTE, repeats=REPEATS, masker_file=MASKER_FILE, max_correlations=MAX_CORRELATIONS, keep_temp=KEEP_TEMP, experiment_b_path=EXPERIMENT,
                 replica_a_path=EXPERIMENT, replica_b_path=EXPERIMENT, poisson_test=POISSONTEST, stranded_analysis=STRANDED_ANALYSIS, proximity=PROXIMITY, 
                 postscript=POSTSCRIPT, showplots=SHOWPLOTS, plot_path=PLOT_PATH, pseudocount=PSEUDOCOUNT, len_norm=LEN_NORM, label1=LABEL1, 
                 label2=LABEL2, binsize=BINSIZE, zscore=ZSCORE, blacklist=BLACKLIST, sdfold=SDFOLD, recalculate=RECALCULATE, counts_file=COUNTS_FILE, 
                 region_mintags=REGION_MINTAGS, bin_step=WINDOW_STEP, tmm_norm=TMM_NORM, n_norm=N_NORM, skip_header=SKIP_HEADER, total_reads_a=TOTAL_READS_A, 
                 total_reads_b=TOTAL_READS_B, total_reads_replica=TOTAL_READS_REPLICA, a_trim=A_TRIM, m_trim=M_TRIM, use_replica_flag=USE_REPLICA, tempdir=TEMPDIR,
                 use_samtools=USESAMTOOLS):

        self.__dict__.update(locals())
        if type(self.tempdir) is not list: #
            self.tempdir = [self.tempdir] 


        self.normalize_factor = 1
        self.set_logger(verbose, debug)
        self.is_sorted = False
        self.temp_experiment = False #Indicates if temp files where created for the experiment
        self.temp_control = False #Indicates if temporary files where created for the control
        self.temp_replica = False 
        self.operations = []
        self.previous_chr = None
        self.open_region = False
        self.safe_reader = utils.SafeReader(self.logger)        
        if not self.discarded_names:
            self.discarded_names = []
     
        #Reusable cluster objects 
        try:
            self.cluster = Cluster(read=self.experiment_format, write=self.output_format, rounding=self.rounding, read_half_open = self.open_experiment, write_half_open = self.open_output, tag_length=self.tag_length, span = self.span, logger=self.logger, cached=cached)
            self.cluster_aux = Cluster(read=self.experiment_format, write=self.output_format, rounding=self.rounding, read_half_open = self.open_experiment, write_half_open = self.open_output, tag_length=self.tag_length, span = self.span, logger=self.logger, cached=cached)
            self.cluster_aux2 = Cluster(read=self.experiment_format, write=self.output_format, rounding=self.rounding, read_half_open = self.open_experiment, write_half_open = self.open_output, tag_length=self.tag_length, span = self.span, logger=self.logger, cached=cached)
            self.control_cluster = Cluster(read=experiment_format, read_half_open = open_experiment, logger=self.logger, cached=cached)

        except ConversionNotSupported:
            self.logger.error('Reading "%s" and writing "%s" is not supported.\n'%(self.experiment_format, self.output_format))
            utils.list_available_formats()

        #duplicates flags
        self.previous_start = 0
        self.previous_end = 0
        self.previous_chr = ''
        self.duplicates_found = 0
        #poisson stuff
        self.first_chr = True
        self._init_poisson()
        self.poisson_results = {}
        for test in POISSON_OPTIONS:
            self.poisson_results[test] = defaultdict()
        
        #self.poisson_results = {'length': defaultdict(), 'height': defaultdict(), 'numtags': defaultdict()}
        self.maxheight_to_pvalue = {}
        self.numtags_to_pvalue = {}
        #Operation flags
        self.do_poisson = False
        self.do_subtract = False
        self.do_heuremove = False
        self.do_split = False
        self.do_trim = False
        self.do_cut = False
        self.do_dupremove = False
        self.sorted_by_pyicos = False
        self.sorted_region_path = ''

    def i_cant_do(self):
        """Quits and exists if exits if the combination of non possible operations"""
        if FILTER in self.operations and POISSON not in self.operations and not self.threshold:
            self.logger.error("Can't do Filter without Poisson or a fixed threshold\n")
            sys.exit(0)
        elif (SUBTRACT or SPLIT or POISSON) in self.operations and (self.output_format not in CLUSTER_FORMATS):
            self.logger.error("Can't get the output as tag format (eland, bed) for Subtract, Split, Poisson filtering please use a clustered format")
            print CLUSTER_FORMATS
            sys.exit(0)
        elif EXTEND in self.operations and self.experiment_format in CLUSTER_FORMATS:
            self.logger.error("Can't extend if the experiment is a clustered format ")
            print CLUSTER_FORMATS
            sys.exit(0)
        elif STRAND_CORRELATION in self.operations and self.experiment_format in CLUSTER_FORMATS:
            self.logger.error("Can't perform strand correlation operation if the experiment is a clustered format ")
            print CLUSTER_FORMATS
            sys.exit(0)
        
        elif ENRICHMENT in self.operations and MODFDR in self.operations:
            print
            self.logger.error("Enrichment and ModFDR operations are not compatible")
            sys.exit(0)

        elif ENRICHMENT in self.operations and POISSON in self.operations:
            print
            self.logger.error("Enrichment and Poisson operations are not compatible")
            sys.exit(0)             
    
        
    def _add_slash_to_path(self, path):
        return utils.add_slash_to_path(path)

    def success_message(self, output_path):
        self.result_log.write('\nSummary of operations\n')
        self.result_log.write('---------------------\n')
        self.logger.info('Operations complete.')
        self.logger.info('')
        self.logger.info('Summary of operations:')
        self.logger.info('----------------------')
        if self.experiment_format != self.output_format and not NOWRITE in self.operations:
            self.logger.info('Convert from %s to %s.'%(self.experiment_format, self.output_format))

        if STRAND_CORRELATION in self.operations:
            self.result_log.write('Strand correlation: %.0f pairs analyzed, estimated a %s nucleotides extension value'%(self.analyzed_pairs, self.frag_size))
            self.logger.info('Strand correlation')

        if EXTEND in self.operations:
            self.result_log.write('Extend to %s\n'%self.frag_size)
            self.logger.info('Extend to %s'%self.frag_size)

        if self.do_subtract:
            self.result_log.write('Subtract\n')
            self.logger.info('Subtract')

        if self.discarded_names:
            self.result_log.write('Discard tags: %s\n'%self.discarded_names)
            self.logger.info('Discard tags: %s'%self.discarded_names)

        if self.do_heuremove:
            self.result_log.write('Heuristic Remove from %s\n'%self.region_path)
            self.logger.info('Heuristic Remove from %s'%self.region_path)

        if self.do_split:
            self.result_log.write('Split\n')
            self.logger.info('Split')

        if DISCARD_ARTIFACTS in self.operations:
            self.result_log.write('Discard Artifacts\n')
            self.logger.info('Discard Artifacts')

        if self.do_poisson:
            self.result_log.write('Poisson\n')
            self.logger.info('Poisson')

        if self.do_cut:
            self.result_log.write('Filter\n')
            self.logger.info('Filter')

        if self.do_dupremove:
            self.result_log.write('Removed %s duplicates\n'%self.duplicates_found)
            self.logger.info('Removed %s duplicates, allowing up to %s'%(self.duplicates_found, self.duplicates))

        self.result_log.write('Date finished: %s'%datetime.now())        
        if not NOWRITE in self.operations:
            self.logger.info('Output at: %s'%(output_path))

    def start_operation_message(self):
        if self.experiment_format != self.output_format and not NOWRITE in self.operations:
            self.logger.info('Converting file %s from %s to %s...'%(self.current_experiment_path, self.experiment_format, self.output_format))
        else:
            self.logger.info('Reading file %s as a %s file...'%(self.current_experiment_path, self.experiment_format))

        if self.current_control_path:
            self.logger.info('Control file:%s'%self.current_control_path)
            if self.do_normalize:
                self.logger.info('The file %s will be normalized to match %s'%(self.current_experiment_path, self.current_control_path))
            if self.do_subtract:
                self.logger.info('The file %s will be subtracted from %s'%(self.current_control_path, self.current_experiment_path))

    def get_normalize_factor(self, experiment, control):
        ret = self.numcells(control, self.control_format)/self.numcells(experiment, self.experiment_format)
        self.result_log.write('Normalization factor: %s\n'%(ret))
        self.logger.info('Normalization factor: %s\n'%(ret))
        return ret
    
    def numcells(self, file_path, file_format):
        """Returns the total number of cells in a file"""
        num = 0.
        cluster = Cluster(read=file_format, logger=self.logger)
        for line in utils.open_file(file_path, file_format, logger=self.logger):
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
        self.do_heuremove = (REMOVE_REGION in self.operations and self.blacklist)
        self.do_poisson = POISSON in self.operations
        self.do_split = SPLIT in self.operations
        self.do_trim = TRIM in self.operations
        self.do_cut = FILTER in self.operations
        self.do_discard = DISCARD_ARTIFACTS in self.operations
        self.do_dupremove = REMOVE_DUPLICATES in self.operations
        self.tag_to_cluster = (not self.experiment_format in CLUSTER_FORMATS and self.output_format in CLUSTER_FORMATS)
        self.use_MA = USE_MA in self.operations
        self.logger.info("")
        self.logger.info('***********************************************')
        self.logger.info('Pyicos running... (PID: %s)'%os.getpid())
        if self.control_path:
            self.process_all_files_paired(self.experiment_path, self.control_path)
        elif self.experiment_b_path:
            self.process_all_files_paired(self.experiment_path, self.experiment_b_path)
        elif self.plot_path:
            self.process_all_files_recursive(self.plot_path, self.output_path)
        elif self.counts_file:
            self.process_all_files_recursive(self.counts_file, self.output_path)
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
            self.logger.info('Calculating normalization factor...')
            self.normalize_factor = self.get_normalize_factor(self.current_experiment_path, self.current_control_path)
            self.cluster.normalize_factor = self.normalize_factor
            self.cluster_aux.normalize_factor = self.normalize_factor
            self.cluster_aux2.normalize_factor = self.normalize_factor

    def _manage_temp_file(self, path):
        """A temporary file that is no longer needed is given, and depending on the value of self.keep_temp its removed or kept"""
        if self.keep_temp:
            self.logger.info("Temp file kept at: %s (remember cleaning them!)"%path)
        else:
            os.remove(path)
            self.logger.info('Temporary file %s removed'%path)

    def _filter_file(self, file_path, temp_name, remove_temp, file_format, file_open):
        """Assumes sorted file. Extend and removal of duplicates go here"""
        if self.tempdir:
            td = self.tempdir[0]
        else:
            td = gettempdir()

        new_file_path = "%s/tempfilter%s_%s"%(td, os.getpid(), temp_name)
        new_file = open(new_file_path, 'w')
        self.logger.info("Filtering %s file..."%temp_name)
        previous_line = ''
        equal_lines = 0
        self.cluster.clear()

        cluster_same = Cluster(read=file_format, write=file_format, rounding=self.rounding, read_half_open = file_open, write_half_open = file_open, tag_length=self.tag_length, span = self.span, logger=self.logger, cached=self.cached)

        for line in utils.open_file(file_path, file_format, logger=self.logger):
            self.cluster.clear()
            self.safe_read_line(self.cluster, line)
            if self.cluster.start == cluster_same.start and self.cluster.end == cluster_same.end and self.cluster.name == cluster_same.name and self.cluster.strand == cluster_same.strand:
                equal_lines+=1
            else:
                equal_lines=0

            cluster_same.clear()
            self.safe_read_line(cluster_same, line)
     
            if self.do_heuremove:
                cluster_same = self.remove_regions(cluster_same)

            if self.duplicates >= equal_lines:
                if EXTEND in self.operations:
                    cluster_same.extend(self.frag_size)

                new_file.write(cluster_same.write_line())
            else:
                self.duplicates_found += 1

            cluster_same.clear() #Need to re-read to allow combination of extend and remove duplicates
            self.safe_read_line(cluster_same, line)

        new_file.flush()
        new_file.close()

        if EXTEND in self.operations:
            self.logger.info("Resorting %s after extension..."%temp_name)
            sorter = utils.BigSort(file_format, file_open, 0, 'fisort%s'%temp_name, logger=self.logger)
            old_path = new_file_path
            sorted_new = sorter.sort(old_path, None, self.get_lambda_func(file_format))
            new_file_path = sorted_new.name
            self._manage_temp_file(old_path)

        if remove_temp:
            self._manage_temp_file(file_path) 

        self.cluster.clear()
        self.cluster_aux.clear()   
        return new_file_path     
    
    def filter_files(self):
        """Filter the files removing the duplicates."""
        #TODO normalize should go here too
        self.current_experiment_path = self._filter_file(self.current_experiment_path, "experiment", self.temp_experiment, self.experiment_format, self.open_experiment)
        if self.experiment_format == BAM: 
            self.experiment_format = SAM
        self.temp_experiment = True
        if self.current_control_path:
            self.current_control_path = self._filter_file(self.current_control_path, "control", self.temp_control, self.control_format, self.open_control)
            self.temp_control = True
            if self.control_format == BAM: 
                self.control_format = SAM

    def get_lambda_func(self, format):
        if self.output_format == SPK:
            if format == ELAND:
                return lambda x:(x.split()[6],x.split()[8],int(x.split()[7]),len(x.split()[1]))
            elif format == SAM or format == BAM:
                return lambda x:(x.split()[2], x.split()[1], int(x.split()[3]), len(x.split()[9]))
            else:
                return lambda x:(x.split()[0], x.split()[5], int(x.split()[1]), int(x.split()[2]))
        else:
            if format == ELAND:
                return lambda x:(x.split()[6], int(x.split()[7]), len(x.split()[1]))
            elif format == SAM or format == BAM:
                return lambda x:(x.split()[2], int(x.split()[3]), len(x.split()[9]))
            else:
                return lambda x:(x.split()[0],int(x.split()[1]),int(x.split()[2]))


    def decide_sort(self, experiment_path, control_path=None):
        """Decide if the files need to be sorted or not."""
        #TODO refractor this, copy pasted code (warning, its not as easy as it seems)
        if self.tag_to_cluster or self.do_subtract or self.do_heuremove or self.do_dupremove or MODFDR in self.operations or ENRICHMENT in self.operations or REMOVE_REGION in self.operations or STRAND_CORRELATION in self.operations or self.frag_size:
            self.experiment_preprocessor = utils.BigSort(self.experiment_format, self.open_experiment, self.frag_size, 'experiment', logger=self.logger)
            self.experiment_b_preprocessor = utils.BigSort(self.experiment_format, self.open_experiment, self.frag_size, 'experiment_b', logger=self.logger)
            self.replica_a_preprocessor = utils.BigSort(self.experiment_format, self.open_experiment, self.frag_size, 'replica_a', logger=self.logger)
            self.control_preprocessor = utils.BigSort(self.control_format, self.open_control, self.frag_size, 'control', logger=self.logger)
                       
            if self.no_sort or self.experiment_format == BAM:
                if self.no_sort: self.logger.warning('Input sort skipped. Results might be wrong.')
                self.current_experiment_path = experiment_path
            elif not self.counts_file:
                self.logger.info('Sorting experiment file...')
                self.is_sorted = True
                sorted_experiment_file = self.experiment_preprocessor.sort(experiment_path, None, self.get_lambda_func(self.experiment_format), self.tempdir)
                self.current_experiment_path = sorted_experiment_file.name
                self.temp_experiment = True

            if self.do_subtract:
                if self.no_sort or self.experiment_format == BAM:
                    if self.no_sort: self.logger.warning('Control sort skipped. Results might be wrong.')
                    self.current_control_path = control_path
                else:
                    self.logger.info('Sorting control file...')
                    sorted_control_file = self.control_preprocessor.sort(control_path, None, self.get_lambda_func(self.control_format), self.tempdir)
                    self.current_control_path = sorted_control_file.name
                    self.temp_control = True
            
            if ENRICHMENT in self.operations:  
                if self.no_sort or self.experiment_format == BAM:
                    if self.no_sort:  self.logger.warning('Experiment_b sort skipped. Results might be wrong.')
                    self.current_control_path = control_path
                elif not self.counts_file:
                    self.logger.info('Sorting experiment_b file...')
                    sorted_control_file = self.experiment_b_preprocessor.sort(control_path, None, self.get_lambda_func(self.experiment_format), self.tempdir)
                    self.current_control_path = sorted_control_file.name
                    self.temp_control = True
            
            if self.replica_a_path:  
                if self.no_sort or self.experiment_format == BAM:
                    if self.no_sort: self.logger.warning('replica_a sort skipped')
                    self.current_replica_a_path = self.replica_a_path
                else:
                    self.logger.info('Sorting replica_a file...')
                    sorted_replica_a_file = self.replica_a_preprocessor.sort(self.replica_a_path, None, self.get_lambda_func(self.experiment_format), self.tempdir)
                    self.current_replica_a_path = sorted_replica_a_file.name
                    self.temp_replica = True

        if self.region_path:
            if self.no_sort or self.experiment_format == BAM:
                if self.no_sort:
                    self.logger.warning('region sort skipped')
                self.sorted_region_path = self.region_path
            else:
                self.logger.info("Sorting region file...")
                self.region_preprocessor = utils.BigSort(self.region_format, self.open_region, None, 'region', logger=self.logger)
                self.sorted_region_file = self.region_preprocessor.sort(self.region_path, None, self.get_lambda_func(BED), self.tempdir)
                self.sorted_region_path = self.sorted_region_file.name



    def operate(self, experiment_path, control_path=None, output_path=None):
        """Operate expects single paths, not directories. Its called from run() several times if the experiment for picos is a directory"""
        #TODO This is where the combination of operations should be done. Operations should be executed in order of inclusion to the list while taking into consideration dependencies and co-runs
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
            self.result_log = open('%s/pyicos_report_%s.txt'%(self._output_dir(), os.path.basename(self.current_output_path)), 'wb')
            self.result_log.write('Pyicos analysis report\n')
            self.result_log.write('----------------------\n\n')
            self.result_log.write('Date run: %s\n'%datetime.now())
            self.result_log.write('Experiment file: %s\n'%self.current_experiment_path)
            if self.current_control_path:
                self.result_log.write('Control file: %s\n'%self.current_control_path)
            self.result_log.write('\n\n')

            self.start_operation_message()
            self.decide_sort(experiment_path, control_path)

            self.estimate_frag_size = self.do_poisson and not self.frag_size
         
            if self.control_path:
                self.control_reader_simple = utils.SortedFileReader(self.current_control_path, self.control_format, rounding=self.rounding, logger=self.logger)
                self.control_reader = utils.read_fetcher(self.current_control_path, self.control_format, rounding=self.rounding, cached=self.cached, logger=self.logger)

            if self.blacklist:
                self.blacklist_reader = utils.read_fetcher(self.blacklist, BED, rounding=self.rounding, cached=self.cached, logger=self.logger)

            if STRAND_CORRELATION in self.operations:
                self.strand_correlation()          
            if self.do_dupremove or (EXTEND in self.operations and (STRAND_CORRELATION in self.operations or self.experiment_format == BAM or self.control_format == BAM)) or self.do_heuremove:
                self.filter_files()

            if NORMALIZE in self.operations:
                self.normalize()


            if self.do_cut: #if we cut, we will round later
                self.cluster.rounding = False
                self.cluster_aux.rounding = False

            if ((not NOWRITE in self.operations) or (NOWRITE in self.operations and POISSON in self.operations)) and not ENRICHMENT in self.operations: 
                self.process_file()
            
            if self.do_poisson: #extract info from the last name and print the thresholds
                self.poisson_analysis(self.previous_chr)
                self.logger.info('\nCluster thresholds for p-value %s:'%self.p_value)
                self.result_log.write('\nCluster thresholds for p-value %s:\n'%self.p_value)
                for name, k in self.poisson_results[self.poisson_test].items():
                    self.logger.info('%s: %s'%(name, k))
                    self.result_log.write('%s: %s\n'%(name, k))

            #Mutually exclusive final operations
            if ENRICHMENT in self.operations:
                self.plot_path = self.enrichment()
                if CALCZSCORE in self.operations:
                    self.plot_path = self.calculate_zscore(self.plot_path)

            elif self.do_cut: 
                self.cut()

            elif MODFDR in self.operations:
                self.modfdr()

            #Plot and summary operations
            if PLOT in self.operations:
                  self.logger.info("Plotting...")
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

            try:
                if self.temp_replica:
                    self._manage_temp_file(self.current_replica_a_path)

            except AttributeError, OSError:
                pass


    def _to_read_conversion(self, experiment, output):
        for line in experiment:
            try:
                self.safe_read_line(self.cluster, line)
                if not self.cluster.is_empty():
                    self.process_cluster(self.cluster, output)
                self.cluster.clear()
            except InsufficientData:
                self.logger.warning('For this type of conversion (%s to %s) you need to specify the tag length with the --tag-length flag'%(self.experiment_format, self.output_format))
                sys.exit(0)


    def _to_cluster_conversion(self, experiment, output):
        #load the first read
        self.logger.debug("_to_cluster_conversion: running clustering...")
        while self.cluster.is_empty():
            self.safe_read_line(self.cluster, experiment.next())

        for line in experiment:
            self.cluster_aux.clear()
            self.safe_read_line(self.cluster_aux, line)
            if not self.cluster_aux.is_empty():
                if self.cluster_aux.touches(self.cluster):
                    self.cluster += self.cluster_aux
                else:
                    if not self.cluster.is_empty():
                        self.process_cluster(self.cluster, output)
                        self.cluster.clear()
                        
                    self.safe_read_line(self.cluster, line)

        if not self.cluster.is_empty(): #Process the last cluster
            self.process_cluster(self.cluster, output)
            self.cluster.clear()

        self.logger.debug("_to_cluster_conversion: Done clustering.")

    def process_file(self):
        self.cluster.clear()
        self.cluster_aux.clear()
        if NOWRITE in self.operations:
            output = None
        else:
            output = open(self.current_output_path, 'wb')

        experiment = utils.open_file(self.current_experiment_path, self.experiment_format, logger=self.logger)
        if self.output_format == WIG or self.output_format == VARIABLE_WIG:
            output.write('track type=wiggle_0\tname="%s"\tvisibility=full\n'%self.label)

        if self.output_format in CLUSTER_FORMATS and not MODFDR in self.operations and not ENRICHMENT in self.operations:
            if not self.experiment_format in CLUSTER_FORMATS:
                self.logger.info('Clustering reads...')
            self._to_cluster_conversion(experiment, output)
        else:
            self._to_read_conversion(experiment, output)

        if not NOWRITE in self.operations:
            output.flush()
            output.close()

    def subtract(self, cluster):
        region = Region(cluster.name, cluster.start, cluster.end, logger=self.logger)
        self.logger.debug("Getting overlapping clusters...")
        self.cluster_aux.clear()
        for c in self.control_reader_simple.overlapping_clusters(region, overlap=EPSILON):
            if self.cluster_aux.is_empty(): #first one
               c.copy_cluster_data(self.cluster_aux)
           
            elif self.cluster_aux.touches(c):
                self.cluster_aux += c

            else:
                region.clusters.append(self.cluster_aux.copy_cluster())
                self.logger.debug("Region clusters length: %s "%len(region.clusters))
                c.copy_cluster_data(self.cluster_aux)

        if self.cluster_aux: #last one
            region.clusters.append(self.cluster_aux.copy_cluster())
            self.logger.debug("Region clusters length: %s Last cluster cache: %s"%(len(region.clusters), region.clusters[-1]._tag_cache))
            self.cluster_aux.clear()

        self.logger.debug("Done adding tags!")
        self.logger.debug("Clustering control (metacluster)...")
        meta = region.get_metacluster()
        del region
        cluster -= meta
        self.logger.debug("OPERATION: Done subtracting (for real)")       
        return cluster

    def remove_regions(self, cluster):
        region = Region(cluster.name, cluster.start, cluster.end)
        if self.blacklist_reader.get_overlaping_clusters(region, overlap=EPSILON):
            cluster.clear() 
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

    def poisson_analysis(self, name=''):
        """
        We do 3 different poisson statistical tests per name for each experiment file:
 
        Cluster analysis:
        This analysis takes as a basic unit the "cluster" profile and performs a poisson taking into account the height
        of the profile. This will help us to better know witch clusters are statistically significant and which are product of chromatine noise

        We do this by calculating the average height of clusters in a given chromosome. Given this mean, we calculate the p-value for each height k the poisson
        function gives us the probability p of one cluster having a height k by chance. 
        With this if we want, for example, to know what is the probability of getting a cluster higher than k = 7, we accumulate the p-values 
        for poisson(k"0..6", mean), and calculate 1 - sum(poisson(k, mean))
        
        Nucleotide analysis:
        This analysis takes the nucleotide as the unit to analize. We give a p-value for each "height"
        of read per nucleotides using an accumulated poisson. With this test we can infer more accurately 
        what nucleotides in the cluster are part of the DNA binding site.


        Number of reads analysis:
        We analize the number of reads of the cluster. Number of reads = sum(xi *yi ) / read_length

        """
        if os.path.isdir(self.output_path):
           out = self.output_path
        else:
           out = os.path.dirname(os.path.abspath(self.output_path))

        #search for the name length in the chr len files
        found = False
        try:
            #print '%s/../chrdesc/%s'%(os.path.dirname(__file__), self.species)
            for line in open('%s/../chrdesc/%s'%(os.path.dirname(__file__), self.species)):
                
                chrom, length = line.split()
                if chrom == name:
                    found = True
                    self.chr_length = int(length)
        except IOError:
            pass #file not found, the warning will be printed
        
        if not found: self.logger.warning("The file containing %s length for assembly %s could not be found, an aproximation will be used"%(name, self.species))
        self.result_log.write('---------------\n')
        self.result_log.write('%s\n'%(name))
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
            
            if name not in self.poisson_results['length'].keys() and p_nucleotide < self.p_value: #if we don't have a height k that is over the p_value yet, write it.
                self.poisson_results["length"][name] = k

            if name not in self.poisson_results['height'].keys() and p_cluster < self.p_value:
                self.poisson_results["height"][name] = k
            
            if name not in self.poisson_results['numtags'].keys() and p_numtags < self.p_value:
                self.poisson_results["numtags"][name] = k

            if k not in self.maxheight_to_pvalue:
                self.maxheight_to_pvalue[k] = {}
            self.maxheight_to_pvalue[k][name] = p_cluster

            if k not in self.numtags_to_pvalue:
                self.numtags_to_pvalue[k] = {}
            self.numtags_to_pvalue[k][name] = p_numtags


            k+=1

    def poisson_retrieve_data(self, cluster):

        if self.estimate_frag_size: #TODO We need to estimate the fragment size if not provided
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

        if cluster.name not in self.discarded_names and not cluster.is_empty():
            if self.previous_chr != cluster.name: #A new name has been detected
                if self.is_sorted:
                    self.logger.info('Processing %s...'%cluster.name)
                    sys.stdout.flush()
                if self.do_poisson and not self.first_chr:
                    self.poisson_analysis(self.previous_chr)
                    self._init_poisson()
                self.previous_chr = cluster.name
                if self.output_format == VARIABLE_WIG:
                    output.write('variableStep\tchrom=%s\tspan=%s\n'%(self.previous_chr, self.span))
                self.first_chr = False

            if self.do_subtract:
                cluster = self.subtract(cluster)
                self.logger.debug("Absolute split and post subtract...")
                for cluster in cluster.absolute_split(threshold=0):
                    self._post_subtract_process_cluster(cluster, output)
            else:
                self._post_subtract_process_cluster(cluster, output)

            #self.logger.debug("Finished process_cluster")

    def _post_subtract_process_cluster(self, cluster, output):
        if self.do_poisson:
            self.poisson_retrieve_data(cluster)

        if not (cluster.is_artifact() and DISCARD_ARTIFACTS in self.operations) and not NOWRITE in self.operations:
            if self.do_trim:
                cluster.trim(self.trim_percentage, self.trim_absolute)

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
        move(os.path.abspath(self.current_output_path), old_output)
        filtered_output = open(self.current_output_path, 'w+')

        if self.poisson_test == 'height':
            msg = "Filtering using cluster height..."
        else:
            msg = "number of reads per cluster..."

        self.logger.info(msg)
        unfiltered_output = open('%s/unfiltered_%s'%(current_directory, os.path.basename(self.current_output_path)), 'w+')
        if self.output_format == WIG or self.output_format == VARIABLE_WIG:
            wig_header = 'track type=wiggle_0\tname="%s"\tvisibility=full\n'%self.label
            filtered_output.write(wig_header)
            unfiltered_output.write(wig_header)
        cut_cluster = Cluster(read=self.output_format, write=self.output_format, rounding=self.rounding, read_half_open = self.open_output, write_half_open = self.open_output, tag_length=self.tag_length, span = self.span, logger=self.logger, normalize_factor = self.normalize_factor)
        self.logger.info('Writing filtered and unfiltered file...')
        for line in open(old_output):
            cut_cluster.clear()
            self.safe_read_line(cut_cluster, line)
            try:
                if self.poisson_test == 'height':
                    cut_cluster.p_value = self.maxheight_to_pvalue[int(round(cut_cluster.max_height()))][cut_cluster.name]
                    
                else:
                    cut_cluster.p_value = self.numtags_to_pvalue[int(round(cut_cluster.area()/self.frag_size))][cut_cluster.name]
            
            except KeyError:
                cut_cluster.p_value = 0 #If the cluster is not in the dictionary, it means its too big, so the p_value will be 0

            try:
                if self.threshold:
                    thres = self.threshold
                else:
                    thres = self.poisson_results[self.poisson_test][cut_cluster.name]
                    
                if cut_cluster.is_significant(thres, self.poisson_test, self.frag_size):
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
            try:
                strand = None
                exome_size = 0
                if self.stranded_analysis:
                    strand = sline[5]
                elif self.region_format == BED12:
                    for size in sline[10].split(","):
                       exome_size += int(size)

                return Region(sline[0], int(sline[1]), int(sline[2]), name2=sline[3], strand=strand, exome_size=exome_size)

            except ValueError:
                pass #discarding header


    def _region_from_dual(self, line):
            try:

                self.cluster_aux.clear()
                self.cluster_aux.read_line(line)
                strand = None
                if self.stranded_analysis:
                    strand = self.cluster_aux.strand
                ret = Region(self.cluster_aux.name, self.cluster_aux.start, self.cluster_aux.end, name2=self.cluster_aux.name2, strand=strand)
                self.cluster_aux.clear()
                return ret                
            except ValueError:
                pass #discarding header




    def __calc_reg_write(self, region_file, count, calculated_region):
        if count > self.region_mintags:
            region_file.write(calculated_region.write()) 


    def calculate_region(self):
        """
        Calculate a region file using the reads present in the both main files to analyze. 
        """
        self.logger.info('Generating regions...')
        self.sorted_region_path = '%s/calcregion_%s.txt'%(self._output_dir(), os.path.basename(self.current_output_path))
        region_file = open(self.sorted_region_path, 'wb')
        dual_reader = utils.DualSortedReader(self.current_experiment_path, self.current_control_path, self.experiment_format, self.logger) 
        if self.stranded_analysis:
            self.calculate_region_stranded(dual_reader, region_file)    
        else:
            self.calculate_region_notstranded(dual_reader, region_file)    

        region_file.flush()

    def __cr_append(self, regions, region):
        regions.append(region)

    def calculate_region_notstranded(self, dual_reader, region_file):
        calculated_region = Region()        

        readcount = 1
        for line in dual_reader:
            if not calculated_region: #first region only
                calculated_region = self._region_from_dual(line)
                calculated_region.end += self.proximity
            else:
                new_region = self._region_from_dual(line)
                new_region.end += self.proximity
                if calculated_region.overlap(new_region):
                    calculated_region.join(new_region)
                    readcount += 1
                else:
                    calculated_region.end -= self.proximity
                    self.__calc_reg_write(region_file, readcount, calculated_region)                 
                    calculated_region = new_region.copy()
                    readcount = 1

        if calculated_region:
            calculated_region.end -= self.proximity
            self.__calc_reg_write(region_file, readcount, calculated_region)                         


    def calculate_region_stranded(self, dual_reader, region_file):
        temp_region_file = open(self.sorted_region_path, 'wb')
        region_plus = Region()
        region_minus = Region()
        regions = []
        numreads_plus = 1
        numreads_minus = 1
        dual_reader = utils.DualSortedReader(self.current_experiment_path, self.current_control_path, self.experiment_format, self.logger)
        for line in dual_reader:
            new_region = self._region_from_dual(line)
            new_region.end += self.proximity
            if not (region_plus and new_region.strand == PLUS_STRAND):
                region_plus = self._region_from_dual(line)
               
            elif not (region_plus and new_region.strand == PLUS_STRAND):
                region_minus = self._region_from_dual(line)

            else:
                if region_plus.overlap(new_region) and region_plus.strand == new_region.strand:
                    region_plus.join(new_region)
                    numreads_plus += 1
                elif region_minus.overlap(new_region) and region_minus.strand == new_region.strand:
                    region_minus.join(new_region)
                    numreads_minus += 1
                else:
                    if new_region.strand == region_plus.strand:
                        region_plus.end -= self.proximity
                        self.__calc_reg_write(region_file, numreads_plus, region_plus)                      
                        region_plus = new_region.copy()  
                        numreads_plus = 1    
                    else:
                        region_minus.end -= self.proximity
                        self.__calc_reg_write(region_file, numreads_minus, region_minus)         
                        region_minus = new_region.copy()  
                        numreads_minus = 1    

        if region_plus:
            region_plus.end -= self.proximity
            regions.append(region_plus)

        if region_minus:
            region_minus.end -= self.proximity
            regions.append(region_minus)

        regions.sort(key=lambda x:(x.name, x.start, x.end, x.strand))
        for region in regions:
            region_file.write(region.write())  
       

    def get_zscore(self, x, mean, sd):    
        if sd > 0:
            return float(x-mean)/sd
        else:
            return 0 #This points are weird anyway 
        

    def plot_enrichment(self, file_path):
        try:
            if self.postscript:
                import matplotlib
                matplotlib.use("PS")

            from matplotlib.pyplot import plot, hist, show, legend, figure, xlabel, ylabel, subplot, axhline, axis
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
                    label_control = '%s(A) VS %s(A)'%(os.path.basename(self.real_experiment_path), os.path.basename(self.replica_a_path))
                else:
                    label_control = 'Background distribution'

            A = []
            A_prime = []
            M = []
            M_significant = []
            A_significant = []
            M_prime = []
            A_medians = []      
            points = []
            minus_points = []
            all_points = []
            figure(figsize=(8,6))
            biggest_A  = -sys.maxint #for drawing
            smallest_A = sys.maxint #for drawing
            biggest_M = 0 #for drawing

            self.logger.info("Loading table...")
            for line in open(file_path):
                sline = line.split()
                try:
                    enrich = dict(zip(enrichment_keys, sline)) 
                    biggest_A = max(biggest_A, float(enrich["A"]))         
                    smallest_A = min(smallest_A, float(enrich["A"]))   
                    biggest_M = max(biggest_M, abs(float(enrich["M"])))          
        
                    biggest_A = max(biggest_A, float(enrich["A_prime"]))         
                    smallest_A = min(smallest_A, float(enrich["A_prime"]))
                    biggest_M = max(biggest_M, abs(float(enrich["M_prime"])))
               
                    positive_point = self.zscore*float(enrich["sd"])+float(enrich["mean"])
                    negative_point = -self.zscore*float(enrich["sd"])+float(enrich["mean"])
                    A_median = float(enrich["A_median"])
                    all_points.append((A_median, positive_point, negative_point))
                    if abs(float(enrich["zscore"])) < self.zscore:
                        M.append(float(enrich["M"]))
                        A.append(float(enrich["A"]))
                    else:
                        M_significant.append(float(enrich["M"]))                
                        A_significant.append(float(enrich["A"]))  
       
                    M_prime.append(float(enrich["M_prime"]))    
                    A_prime.append(float(enrich["A_prime"]))
                except ValueError:
                    pass #to skip the header
            all_points.sort(key= lambda x:x[0])
            
            for t in all_points:
                (A_medians.append(t[0]), points.append(t[1]), minus_points.append(t[2])) 
                    
            if points:
                margin = 1.1
                A_medians.append(biggest_A*margin)
                points.append(points[-1])
                minus_points.append(minus_points[-1])
                A_medians.insert(0, smallest_A)
                points.insert(0, points[0])
                minus_points.insert(0, minus_points[0])
                self.logger.info("Plotting points...")
                subplot(211)
                xlabel('A')
                ylabel('M')

                axis([smallest_A*margin, biggest_A*margin, -biggest_M*margin, biggest_M*margin])
                plot(A_prime, M_prime, '.', label=label_control, color = '#666666')
                plot(A_medians, points, 'r--', label="z-score %s"%self.zscore)            
                plot(A_medians, minus_points,  'r--')            
                axhline(0, linestyle='--', color="grey", alpha=0.75)  
                legend(bbox_to_anchor=(0., 1.01, 1., .101), loc=3, ncol=1, mode="expand", borderaxespad=0.)

                subplot(212)
                axis([smallest_A*margin, biggest_A*margin, -biggest_M*margin, biggest_M*margin])
                plot(A, M, 'k.', label=label_main)
                plot(A_significant, M_significant, 'r.', label="%s (significant)"%label_main)
                plot(A_medians, points, 'r--', label="z-score %s"%self.zscore)            
                plot(A_medians, minus_points, 'r--')
                axhline(0, linestyle='--', color="grey", alpha=0.75)
                xlabel('A')
                ylabel('M')

                legend(bbox_to_anchor=(0., 1.01, 1., .101), loc=3, ncol=1, mode="expand", borderaxespad=0.)
                self._save_figure("enrichment_MA")
            else:
                self.logger.warning("Nothing to plot.")
        except ImportError:
            if self.debug:
                raise
            self.logger.warning('Pyicos can not find an installation of matplotlib, so no plot will be drawn. If you want to get a plot with the correlation values, install the matplotlib library.')

    def __enrichment_M(self, signal_a, signal_b):
        return math.log(float(signal_a)/float(signal_b), 2)

    def __enrichment_A(self, signal_a, signal_b):
        return (math.log(float(signal_a), 2)+math.log(float(signal_b), 2))/2    
       
         

    def _calculate_MA(self, region_path, read_counts, factor = 1, replica_factor = 1, file_a_reader=None, file_b_reader=None, replica_a_reader=None):
        tags_a = []
        tags_b = []
        numreads_background_1 = 0
        numreads_background_2 = 0
        total_reads_background_1 = 0
        total_reads_background_2 = 0
        self.logger.debug("Inside _calculate_MA")
        self.regions_analyzed_count = 0
        enrichment_result = [] #This will hold the name, start and end of the region, plus the A, M, 'A and 'M
        out_file = open(self.current_output_path, 'wb')
        for region_line in open(region_path):
            sline = region_line.split()

            region_of_interest = self._region_from_sline(sline)            
            if region_of_interest:
                region_a = None
                replica_a = None
                replica_tags = None
                signal_a = -1
                signal_b = -1
                signal_background_1 = -1
                signal_background_2 = -1
                swap1 = Region()
                swap2 = Region()
                if read_counts:
                    signal_a = float(sline[6])
                    signal_b = float(sline[7])*factor
                    signal_background_1 = float(sline[8])
                    signal_background_2 = float(sline[9])*replica_factor
                else:
                    tags_a = file_a_reader.get_overlaping_clusters(region_of_interest, overlap=EPSILON)
                    if len(tags_a) > 100: self.logger.debug("Getting tags...")  
                    tags_b = file_b_reader.get_overlaping_clusters(region_of_interest, overlap=EPSILON)
                    if len(tags_a) > 100: self.logger.debug("tags_a: %s tags_b: %s"%(len(tags_a), len(tags_b)))
                    if self.use_replica:            
                        replica_tags = replica_a_reader.get_overlaping_clusters(region_of_interest, overlap=EPSILON)

                    #if we are using pseudocounts, use the union, use the intersection otherwise
                    if (self.pseudocount and (tags_a or tags_b)) or (not self.pseudocount and tags_a and tags_b): 
                        region_a = region_of_interest.copy()
                        region_b = region_of_interest.copy()
                        region_a.add_tags(tags_a, clusterize=False)
                        region_b.add_tags(tags_b, clusterize=False)
                        signal_a = region_a.normalized_counts(self.len_norm, self.n_norm, self.total_regions, self.pseudocount, factor, self.total_reads_a)
                        signal_b = region_b.normalized_counts(self.len_norm, self.n_norm, self.total_regions, self.pseudocount, factor, self.total_reads_b)
                        self.already_norm = True
  
                
                if not self.counts_file:
                    if (self.pseudocount and (tags_a or tags_b)) or (not self.pseudocount and tags_a and tags_b): 
                        if self.use_replica:
                            replica_a = region_of_interest.copy()
                            
                            replica_a.add_tags(replica_tags)
                            numreads_background_1 = len(tags_a)
                            numreads_background_2 = len(replica_tags)
                            total_reads_background_1 = self.total_reads_a
                            total_reads_background_2 = self.total_reads_replica
                            signal_background_1 = signal_a
                            signal_background_2 = replica_a.normalized_counts(self.len_norm, self.n_norm, self.total_regions, self.pseudocount, replica_factor, self.total_reads_replica)

                        elif region_a:
                            swap1, swap2 = region_a.swap(region_b, factor)
                            if len(tags_a) > 100: self.logger.debug("Swap1 Region: %s len: %s "%(swap1, len(swap1.tags)))
                            numreads_background_1 = len(swap1.tags)
                            numreads_background_2 = len(swap2.tags)
                            total_reads_background_1 = total_reads_background_2 = self.average_total_reads
                            signal_background_1 = swap1.normalized_counts(self.len_norm, self.n_norm, self.total_regions, self.pseudocount, replica_factor, self.average_total_reads)
                            signal_background_2 = swap2.normalized_counts(self.len_norm, self.n_norm, self.total_regions, self.pseudocount, replica_factor, self.average_total_reads)

                #if there is no data in the replica or in the swap and we are not using pseudocounts, dont write the data 
                if signal_a > 0 and signal_b > 0 and signal_background_1 > 0 and signal_background_2 > 0 or self.use_MA:
                    if self.use_MA and not self.already_norm:
                        A = float(sline[10])
                        M = float(sline[11])
                        A_prime = float(sline[16])
                        M_prime = float(sline[17])
                    else:
                        if not self.already_norm: #TODO refractor
                            if self.len_norm: #read per kilobase in region
                                signal_a = 1e3*(float(signal_a)/len(region_of_interest))
                                signal_b = 1e3*(float(signal_b)/len(region_of_interest))
                                signal_background_1 = 1e3*(float(signal_background_1)/len(region_of_interest))
                                signal_background_2 = 1e3*(float(signal_background_2)/len(region_of_interest))
  
                            if self.n_norm: #per million reads in the sample
                                signal_a = 1e6*(float(signal_a)/self.total_reads_a)
                                signal_b = 1e6*(float(signal_b)/self.total_reads_b)
                                if self.use_replica:
                                    signal_background_1 = signal_a
                                    signal_background_2 = 1e6*(float(signal_background_2)/self.total_reads_replica)
                                else:
                                    signal_background_1 = 1e6*(float(signal_background_1)/self.average_total_reads)
                                    signal_background_2 = 1e6*(float(signal_background_2)/self.average_total_reads)                                                                    
                                
                        A = self.__enrichment_A(signal_a, signal_b)
                        M = self.__enrichment_M(signal_a, signal_b)
                        A_prime = self.__enrichment_A(signal_background_1, signal_background_2)
                        M_prime = self.__enrichment_M(signal_background_1, signal_background_2)

                    out_file.write("%s\n"%("\t".join([region_of_interest.write().rstrip("\n"), str(signal_a), str(signal_b), str(signal_background_1), str(signal_background_2), str(A), str(M), str(self.total_reads_a), str(self.total_reads_b), str(len(tags_a)), str(len(tags_b)),  str(A_prime), str(M_prime), str(total_reads_background_1), str(total_reads_background_2), str(numreads_background_1), str(numreads_background_2)])))
                    self.regions_analyzed_count += 1

        out_file.flush()
        out_file.close()
        self.logger.debug("LEAVING _calculate_MA")
        return out_file.name


    def _calculate_total_lengths(self):
        msg = "Calculating enrichment in regions"
        if self.counts_file: 
            self.sorted_region_path = self.counts_file
            if (not self.total_reads_a or not self.total_reads_b or (not self.total_reads_replica and self.use_replica)) and not self.use_MA:
                self.logger.info("... counting from counts file...")
                self.total_reads_a = 0
                self.total_reads_b = 0
                if self.total_reads_replica:
                    self.total_reads_replica = 0
                else:
                    self.total_reads_replica = 1
                for line in open(self.counts_file):
                    try:
                        enrich = dict(zip(enrichment_keys, line.split()))
                        self.total_reads_a += float(enrich["signal_a"])
                        self.total_reads_b += float(enrich["signal_b"])
                        if self.use_replica:
                            self.total_reads_replica += float(enrich["signal_prime_2"])
                    except ValueError:
                        self.logger.debug("(Counting) skip header...")


        else:
            self.logger.info("... counting number of lines in files...")
            if not self.total_reads_a:
                if self.experiment_format == BAM:
                    self.total_reads_a = bam.size(self.current_experiment_path)
                else:
                    self.total_reads_a = sum(1 for line in utils.open_file(self.current_experiment_path, self.experiment_format, logger=self.logger))
            if not self.total_reads_b:
                if self.experiment_format == BAM:
                    self.total_reads_b = bam.size(self.current_control_path)
                else:
                    self.total_reads_b = sum(1 for line in utils.open_file(self.current_control_path, self.control_format, logger=self.logger))
            if self.use_replica and not self.total_reads_replica:
                if self.experiment_format  == BAM:
                    self.total_reads_replica = bam.size(self.replica_a_path)
                else:
                    self.total_reads_replica = sum(1 for line in utils.open_file(self.replica_a_path, self.experiment_format, logger=self.logger))

            if self.use_replica:
                msg = "%s using replicas..."%msg
            else:
                msg = "%s using swap..."%msg

            self.logger.info(msg)

        self.average_total_reads = (self.total_reads_a+self.total_reads_b)/2        


    def enrichment(self):
        file_a_reader = file_b_reader = replica_a_reader = None
        self.use_replica = (bool(self.replica_a_path) or (bool(self.counts_file) and self.use_replica_flag))
        self.logger.debug("Use replica: %s"%self.use_replica)
        
        self._calculate_total_lengths()
        if not self.counts_file:
            file_a_reader = utils.read_fetcher(self.current_experiment_path, self.experiment_format, cached=self.cached, logger=self.logger, use_samtools=self.use_samtools)
            file_b_reader = utils.read_fetcher(self.current_control_path, self.experiment_format, cached=self.cached, logger=self.logger, use_samtools=self.use_samtools)
            if self.use_replica:
                replica_a_reader = utils.read_fetcher(self.current_replica_a_path, self.experiment_format, cached=self.cached, logger=self.logger, use_samtools=self.use_samtools)

            if self.sorted_region_path:
                self.logger.info('Using region file %s (%s)'%(self.region_path, self.region_format))
            else:
                self.calculate_region() #create region file semi automatically

            self.total_regions = sum(1 for line in open(self.sorted_region_path))

        self.logger.info("... analyzing regions, calculating normalized counts, A / M and replica or swap...")

        self.already_norm = False
        out_path = self._calculate_MA(self.sorted_region_path, bool(self.counts_file), 1, 1, file_a_reader, file_b_reader, replica_a_reader)
        self.already_norm = True
        self.logger.debug("Already normalized: %s"%self.already_norm)
        if self.tmm_norm:
            self.logger.info("TMM Normalizing...")
            tmm_factor = self.tmm_factor(out_path, self.regions_analyzed_count, False)
            replica_tmm_factor = 1
            if self.use_replica:
                replica_tmm_factor = self.tmm_factor(out_path, self.regions_analyzed_count, True)
            #move output file to old output
            #use as input
            old_output = '%s/notnormalized_%s'%(self._current_directory(), os.path.basename(self.current_output_path))
            move(os.path.abspath(self.current_output_path), old_output)
            out_path = self._calculate_MA(old_output, True, tmm_factor, replica_tmm_factor, True) #recalculate with the new factor, using the counts again

        self.logger.info("%s regions analyzed."%self.regions_analyzed_count)
        self.logger.info("Enrichment result saved to %s"%self.current_output_path)
        return out_path



    def _sub_tmm(self, counts_a, counts_b, reads_a, reads_b):
        return (counts_a-reads_a)/(counts_a*reads_a) + (counts_b-reads_b)/(counts_b*reads_b)            


    def tmm_factor(self, file_counts, total_regions, replica):
        if replica:
            signal_1 = "signal_prime_1"
            signal_2 = "signal_prime_2"
            M = "M_prime"
            reads_2 = self.total_reads_replica            
        else:
            signal_1 = "signal_a"
            signal_2 = "signal_b"
            M = "M"
            reads_2 = self.total_reads_b

        values_list = []
        #read the file inside the values_list
        for line in open(file_counts):
            sline = line.split()
            values_list.append(dict(zip(enrichment_keys, sline)))

        a_trim_number = int(round(total_regions*self.a_trim))
        #discard the bad A
        self.logger.debug("Removing the worst A (%s regions, %s percent)"%(a_trim_number, self.a_trim*100))

        values_list.sort(key=lambda x:float(x["A"])) #sort by A
        for i in range (0, a_trim_number):
            values_list.pop(0)

        values_list.sort(key=lambda x:float(x[M])) #sort by M              
        m_trim_number = int(round(total_regions*(self.m_trim/2))) #this number is half the value of the flag, because we will trim half below, and half over 

        #remove on the left
        for i in range(0, m_trim_number):
            values_list.pop(0)
        #remove on the right
        for i in range(0, m_trim_number):
            values_list.pop(-1)

        #now calculate the normalization factor
        arriba = 0
        abajo = 0
        for value in values_list:
            w = self._sub_tmm(float(value[signal_1]), float(value[signal_2]), self.total_reads_a, reads_2)                  
            arriba += w*float(value[M])
            abajo += w
        try:
            factor = 2**(arriba/abajo)
        except ZeroDivisionError:
            self.logger.warning("Division by zero, TMM factor could not be calculated.")
            factor = 1

        if replica:    
            self.logger.info("Replica TMM Normalization Factor: %s"%factor)
        else:
            self.logger.info("TMM Normalization Factor: %s"%factor)  
          
        return factor

    def __load_enrichment_result(self, values_path):
        ret = []
        for line in open(values_path):
            sline = line.split()
            try:
                float(sline[1])
                ret.append(dict(zip(enrichment_keys, sline)))
            except ValueError:
                pass

        return ret        
    






    def calculate_zscore(self, values_path): 
        num_regions = sum(1 for line in open(values_path))
        bin_size = int(self.binsize*num_regions)
        if bin_size < 50:
            self.logger.warning("The bin size results in a sliding window smaller than 50, adjusting window to 50 in order to get statistically meaningful results.")
            bin_size = 50

        bin_step = max(1, int(round(self.bin_step*bin_size)))
        self.logger.info("Enrichment window calculation using a sliding window size of %s, sliding with a step of %s"%(bin_size, bin_step))
        self.logger.info("... calculating zscore...")
        enrichment_result = self.__load_enrichment_result(values_path)
        enrichment_result.sort(key= lambda x:(float(x["A_prime"])))
        self.logger.debug("Number of loaded counts: %s"%len(enrichment_result))        
        self.points = []
        #get the standard deviations
        for i in range(0, num_regions-bin_size+bin_step, bin_step):
            #get the slice
            if i+bin_size < num_regions:
                result_chunk = enrichment_result[i:i+bin_size] 
            else:
                result_chunk = enrichment_result[i:]  #last chunk

            #retrieve the values
            mean_acum = 0
            a_acum = 0
            Ms_replica = []
            for entry in result_chunk:
                mean_acum += float(entry["M_prime"])
                a_acum += float(entry["A_prime"])
                Ms_replica.append(float(entry["M_prime"]))

            #add them to the points of mean and sd
            mean = mean_acum/len(result_chunk)
            sd = math.sqrt((sum((x - mean)**2 for x in Ms_replica))/len(Ms_replica))
            #the A median
            A_median = a_acum / len(result_chunk)

            self.points.append([A_median, mean, sd]) #The A asigned to the window, the mean and the standard deviation  
            #self.logger.debug("Window of %s length, with A median: %s mean: %s sd: %s"%(len(result_chunk), self.points[-1][0], self.points[-1][1], self.points[-1][2], len(self.points)))  

        #update z scores
        for entry in enrichment_result:
            entry["A_median"] = 0
            entry["mean"] = 0
            entry["sd"] = 0
            entry["zscore"] = 0
            closest_A = sys.maxint
            sd_position = 0
            for i in range(0, len(self.points)):
                new_A = self.points[i][0]
                if new_A != closest_A: #skip repeated points
                    if abs(closest_A - float(entry["A"])) >= abs(new_A - float(entry["A"])):
                        closest_A = new_A
                        sd_position = i
                    else:
                        break #already found, no need to go further since the points are ordered
                        
            entry["A_median"] = closest_A
            if self.points: #only calculate if there where enough windows, otherwise give a warning
                self.__sub_zscore(entry, self.points[sd_position]) #the value is the biggest, use the last break

        
        if not self.points: 
            self.logger.warning("Insufficient number of regions analyzed (%s), z-score values could not be calculated"%num_regions)

        enrichment_result.sort(key=lambda x:(x["name"], int(x["start"]), int(x["end"])))
        old_file_path = '%s/before_zscore_%s'%(self._current_directory(), os.path.basename(values_path)) #create path for the outdated file
        move(os.path.abspath(values_path), old_file_path) #move the file
        new_file = file(values_path, 'w') #open a new file in the now empty file space  
        if not self.skip_header:
            new_file.write('\t'.join(enrichment_keys))
            new_file.write('\n') 

        for entry in enrichment_result:    
            new_file.write("\t".join(str(entry[key]) for key in enrichment_keys)+"\n")

        self._manage_temp_file(old_file_path)
        return values_path


    def __sub_zscore(self, entry, point):
        entry["mean"] = str(point[1])
        entry["sd"] = str(point[2])
        entry["zscore"] = str(self.get_zscore(float(entry["M"]), float(entry["mean"]), self.sdfold*float(entry["sd"])))


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
            import warnings 
            with warnings.catch_warnings(): #To make the RuntimeWarning go away
                warnings.simplefilter("ignore")
                savefig(figure_path)
            self.logger.info("%s figure saved to %s"%(figure_name, figure_path))
        clf()
        

    def modfdr(self):
        self.logger.info("Running modfdr filter with %s p-value threshold and %s repeats..."%(self.p_value, self.repeats))
        old_output = '%s/before_modfdr_%s'%(self._current_directory(), os.path.basename(self.current_output_path))
        move(os.path.abspath(self.current_output_path), old_output)
        cluster_reader = utils.read_fetcher(old_output, self.output_format, cached=self.cached)
        #if self.masker_file:
        #    masker_reader = utils.read_fetcher(self.masker_file, BED, cached=self.cached)
        filtered_output = open(self.current_output_path, 'w+')
        unfiltered_output = open('%s/unfiltered_%s'%(self._current_directory(), os.path.basename(self.current_output_path)), 'w+')
        for region_line in open(self.sorted_region_path):
            region = self._region_from_sline(region_line.split())
            region.add_tags(cluster_reader.get_overlaping_clusters(region, overlap=EPSILON), True)
            classified_clusters = region.get_FDR_clusters(self.repeats)
            for cluster in classified_clusters:
                unfiltered_output.write(cluster.write_line())
                if cluster.p_value < self.p_value:
                    filtered_output.write(cluster.write_line())

        self._manage_temp_file(old_output)


    def strand_correlation(self):
        self.logger.info("Strand correlation analysis...")
        self.delta_results = dict()
        self.best_delta = -1
        positive_cluster = Cluster(cached=True)
        negative_cluster = Cluster(cached=True)
        positive_cluster_cache = [] #we are trying to hold to the previous cluster
        self.analyzed_pairs = 0.
        acum_length = 0.
        num_analyzed = 0
        for line in utils.open_file(self.current_experiment_path, self.experiment_format, logger=self.logger):
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
                while distance > self.max_delta or positive_cluster.name < negative_cluster.name: # if the negative clusters are too far behind, empty the positive cluster
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
            self.logger.info("Average analyzed length %s"%average_len)
            self.result_log.write("Average analyzed length:%s\n"%average_len)
            for delta in range(self.min_delta, self.max_delta, self.delta_step):
                if delta in self.delta_results:
                    self.delta_results[delta]=self.delta_results[delta]/self.analyzed_pairs
                    data.append(self.delta_results[delta])
                    if self.delta_results[delta] > max_corr:
                        max_delta = delta
                        max_corr = self.delta_results[delta]
                    self.logger.debug('Delta %s:%s'%(delta, self.delta_results[delta]))
        self.logger.info('Correlation test RESULT: You should extend this dataset to %s nucleotides'%(max_delta+average_len))
        self.result_log.write('Correlation test RESULT: You should extend this dataset to %s nucleotides\n'%(max_delta+average_len))
        self.frag_size = int(round(max_delta+average_len))


        if not data:
            self.logger.warning('Not enough data to plot the correlation graph. Lower the threshold of the --height-filter flag')
        else: 
            try:
                if self.postscript:
                    import matplotlib
                    matplotlib.use("PS")

                import matplotlib.pyplot as plt
                plt.plot(range(self.min_delta, self.max_delta), data)
                plt.plot()
                plt.xlabel("Shift length between + and - clusters")
                plt.ylabel("Pearson correlation coefficient")
                self._save_figure("correlation")

            except ImportError:
                if self.debug:
                    raise
                self.logger.warning('Pyicos can not find an installation of matplotlib, so no plot will be drawn for the strand correlation. If you want to get a plot with the correlation values, install the matplotlib library (version >1.0.1).')


    def _correlate_clusters(self, positive_cluster, negative_cluster):
        distance = negative_cluster.end-positive_cluster.start
        if (distance < self.max_delta and distance > self.min_delta) and (positive_cluster.name == negative_cluster.name):
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



