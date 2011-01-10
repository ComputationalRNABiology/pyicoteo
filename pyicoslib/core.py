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

import re
import sys
import random
import math
from collections import defaultdict
from defaults import *

debug = False
verbose = False


class Empty:
    def __init__(self):
        x = 1

##################################
#  EXCEPTIONS                    #
##################################
class InvalidLine(Exception):
    pass

class NotClusteredException(Exception):
    pass

class ConversionNotSupported(Exception):
    pass

class InsufficientData(Exception):
    pass

class DifferentChromosome(Exception):
    pass


class AbstractCore:
    """
    Do not cast, abstract class
    """
    def intersects(self, other):
        """Returns true if a Cluster/Region intersects with another Cluster/Region"""
        return ((other.start >= self.start and other.start <= self.end)
                 or (other.end >= self.start and other.end < self.end)
                 or (other.start <= self.start and other.end >= self.end)) and (self.chromosome == other.chromosome)

    def overlap(self, other):
        """Returns the percentage of overlap of the self cluster with another cluster, from 0 to 1"""
        if self.chromosome != other.chromosome or not self.intersects(other): #different chromosome or no intersection, no overlap
            return 0
        if other.start > self.start and other.end >= self.end:
            #|--------------------| self
            #           |--------------------| other
            #OR
            #|--------------------| self
            #           |---------| other
            return float(self.end-other.start+1)/len(self)

        elif other.end < self.end and other.start <= self.start:
            #      |-------------------| self
            #|---------------| other
            #OR
            #|---------------------| self
            #|---------------| other
            return float(other.end-self.start+1)/len(self)
        elif self.start < other.start and self.end > other.end:
            #|------------------------| self
            #    |------------| other
            return float(len(other))/len(self)
        else:
            #    |------------| self
            #|------------------------| other
            #OR
            #|------------------------| self
            #|------------------------| other
            #OR
            return 1

    def is_contiguous(self, other):
        """Returns true if a Cluster read is contiguous to another one. """
        return (other.start == self.end+1 or other.end+1 == self.start) and self.chromosome == other.chromosome

    

###################################
#   READERS                       #
###################################
class ReaderFactory:
    def create_reader(self, format, half_open=False, cached=True):
        if format == BED:
            return BedReader(format, half_open, cached)
        elif format == PK or format == SPK:
            return PkReader(format, half_open, cached)
        elif format == WIG:
            return WigReader(format, half_open, cached)
        elif format == ELAND:
            return ElandReader(format, half_open, cached)
        elif format == SAM:
            return SamReader(format, half_open, cached)
        else:
            
            raise ConversionNotSupported

class Reader:

    def __init__(self, format, half_open, cached):
        if half_open:
            self.correction = 1
        else:
            self.correction = 0

        self.format = format
        self.half_open = half_open
        self.cached = cached


    def read_line(self):
        raise NotImplementedError("You're using the abstract base class 'Reader', use a specific class instead")

    def _add_line_to_cluster(self, line, cluster):
        cluster_aux = Cluster()
        cluster_aux.reader = cluster.reader
        cluster_aux.normalize_factor = cluster.normalize_factor
        cluster_aux.read_line(line)
        result = cluster + cluster_aux
        cluster.start = result.start
        cluster.end = result.end
        cluster._levels = result._levels
        #if the strands are different, the cluster becomes strandless
        if cluster.strand != cluster_aux.strand:
            cluster.strand == '.'


    def quality_filter(self, line):
        """checks if the line passes the quality conditions"""
        return True

class BedReader(Reader):

    def _add_name(self, cluster, line):
        if len(line) > 3:
            cluster.name = line[3]
        else:
            cluster.name = 'pyicos'

    def _add_score(self, cluster, line):
        if len(line) > 4:
            cluster.score = line[4] #this should be casted to int, but since we dont use it better that its not, lot of people mix this label and the 'name' one
        else:
            cluster.score = 300

    def _add_strand(self, cluster, line):
        if len(line) > 5:
            cluster.strand = line[5]
        else:
            cluster.strand = '.'

    def read_line(self, cluster, line):
        cluster.read_count += 1
        try:
            if self.cached:
                line = line.split()
                new_start = int(line[1])+self.correction
                cluster.chromosome = line[0]
                if cluster.is_empty():
                    cluster._tag_cache.append([new_start, int(line[2])])
                    cluster.start = new_start
                    cluster.end = int(line[2])
                    cluster.tag_length = len(cluster) #if the cluster is empty, its the perfect moment to read the tag length of the first read. We assume that all reads are the same, inside the cluster, so in the case that they vary (not that common), the randomness of choosing this should cancel the bias.  
                    self._add_strand(cluster, line)
                    
                else:
                    cluster._tag_cache.append([new_start, int(line[2])])
                    cluster.start = min(cluster.start, new_start)
                    cluster.end = max(cluster.end, int(line[2]))
                    if len(line) > 5:
                        if line[5] != cluster.strand:
                            cluster.strand = '.'
                    else:
                        cluster.strand = '.'
                
                self._add_name(cluster, line)
                self._add_score(cluster, line)
            else:
                if not cluster.is_empty():
                    self._add_line_to_cluster(line, cluster)
                else:
                    line = line.split()
                    cluster.chromosome = line[0]
                    cluster.start = int(line[1])+self.correction
                    cluster.end = int(line[2])
                    self._add_name(cluster, line)
                    self._add_score(cluster, line)
                    self._add_strand(cluster, line)
                    cluster.tag_length = len(cluster)
                    cluster.add_level(cluster.end-cluster.start+1, cluster.normalize_factor)

        except (ValueError, IndexError):
            raise InvalidLine

class SamReader(Reader):

    def _get_strand(self):
        strand = "+"
        if (self.sam_flag & (0x10)):	# minus strand if true.
            strand = "-"
        return strand

    def read_line(self, cluster, line):
        try:
            line = line.split()
            self.sam_flag = int(line[1])
            if (not (self.sam_flag & 0x0004)):
                new_start = int(line[3])+self.correction
                new_end = new_start+len(line[9])
                cluster.chromosome = line[2]
                new_strand = self._get_strand()
                if cluster.is_empty():
                    cluster._tag_cache.append([new_start, new_start+len(line[9])])
                    cluster.start = new_start
                    cluster.end = new_end
                    cluster.name = line[0]
                    cluster.sequence = line[9]
                    cluster.strand = new_strand
                    cluster.tag_length = len(cluster)
                    
                else:
                    cluster._tag_cache.append([new_start, new_start+len(line[9])])
                    cluster.start = min(cluster.start, new_start)
                    cluster.end = max(cluster.end, new_end)
                    cluster.sequence = ''
                    cluster.name = ''
                    if cluster.strand != new_strand:
                        cluster.strand == '.'

        except (ValueError, IndexError):
            raise InvalidLine

    def quality_filter(self, line):
        sline = line.split()
        try:
            mapped = (int(sline[1]))
        except:
            return False

        return not (mapped & 0x0004)


class WigReader(Reader):
    def read_line(self, cluster, line):
        cluster.read_count = None
        try:
            line = line.split()
            if cluster.is_empty():
                cluster.chromosome = line[0]
                cluster.start = int(line[1])+self.correction

            cluster.add_level(int(line[2])-int(line[1])-self.correction+1, float(line[3])*cluster.normalize_factor)
            cluster._recalculate_end()
        except (ValueError, IndexError):
            raise InvalidLine

class ElandReader(Reader):
    eland_filter = re.compile(r'chr\w{1,2}.fa')

    def read_line(self, cluster, line):
        if line is not None and line != '\n' and self.quality_filter(line):
            try:
                cluster.read_count += 1
                if cluster.is_empty():
                    line = line.split()
                    cluster.name = line[0]
                    length = len(line[1])
                    cluster.sequence = line[1]
                    cluster.chromosome = line[6].rstrip('.fa')
                    cluster.start = int(line[7])+self.correction
                    cluster.end = length+cluster.start-1
                    cluster._levels.append([length+self.correction, cluster.normalize_factor])
                    cluster.tag_length = len(cluster)
                    if line[8] is 'F':
                        cluster.strand = '+'
                    else:
                        cluster.strand = '-'
                else:
                    self._add_line_to_cluster(line, cluster)

            except (ValueError, IndexError):
                raise InvalidLine

    def quality_filter(self, line):
        return self.eland_filter.search(line)

class PkReader(Reader):
    def read_line(self, cluster, line):
        try:
            cluster.read_count = None
            if line is not None and line != '\n':
                if not cluster.is_empty():
                    self._add_line_to_cluster(line, cluster)
                else:
                    line = line.split()
                    if line[0] == 'track':
                        raise InvalidLine
                    cluster.chromosome = line[0]
                    cluster.start = int(line[1])+self.correction
                    for item in line[3].split('|'):
                        temp = item.split(':')
                        cluster.add_level(int(temp[0]), float(temp[1])*cluster.normalize_factor) #normalizing in execution time

                    if len(line) > 5:
                        cluster.strand = line[5]
                    else:
                        cluster.strand = '.'

                    if len(line) > 8:
                        cluster.p_value = float(line[8])

                    cluster._recalculate_end()

        except (ValueError, IndexError):
            raise InvalidLine





##################################
#       WRITERS                  #
##################################
class WriterFactory:
    def create_writer(self, format, half_open, span=0):
        if format == ELAND:
            return ElandWriter(format, half_open, span)
        if format == BED:
            return BedWriter(format, half_open, span)
        elif format == WIG:
            return WigWriter(format, half_open, span)
        elif format == VARIABLE_WIG:
            return VariableWigWriter(format, half_open, span)
        elif format == PK or format == SPK:
            return PkWriter(format, half_open, span)
        elif format == SAM:
            return SamWriter(format, half_open, span)
        else:
            raise ConversionNotSupported

class Writer:
    def __init__(self, format, half_open, span=0):
        if half_open:
            self.correction = -1
        else:
            self.correction = 0
        self.half_open = half_open
        self.span = span
        self.format = format

    def write_line(self):
        raise NotImplementedError("You're using the abstract base class 'Writer', use a specific class instead")

class ElandWriter(Writer):
    def write_line(self, cluster):
        if cluster.is_empty():
            return ''
        else:
            if cluster.strand is '-':
                cluster.strand = 'R'
            else:
                cluster.strand = 'F'
            if cluster.sequence is None:
                seq = 'A'*(cluster.end-cluster.start+1+self.correction)
            else:
                seq = cluster.sequence
            return '%s\t%s\tU0\t1\t0\t0\t%s.fa\t%s\t%s\t..\t26A\n'%(cluster.name, seq, cluster.chromosome, cluster.start+self.correction, cluster.strand)

class BedWriter(Writer):
    def write_line(self, cluster):
        bed_blueprint = '%s\t%s\t%s\t%s\t%s\t%s\n'
        if cluster.is_empty():
            return ''
        else:
            if not cluster.is_singleton():
                lines = ''
                split_tags = cluster.absolute_split(0)
                for tag in split_tags:
                    tags = cluster.decompose()
                    for tag in tags:
                        lines = '%s%s'%(lines, bed_blueprint%(tag.chromosome, tag.start+self.correction, tag.end, tag.name, tag.score, tag.strand))
                return lines
            else:
                return bed_blueprint%(cluster.chromosome, cluster.start+self.correction, cluster.end, cluster.name, cluster.score,  cluster.strand)

class SamWriter(Writer):
    def write_line(self, cluster):
        sam_blueprint = '%s\t%s\t%s\t%s\t255\t%sM0S\t=\t0\t0\t%s\t%s\n'
        if cluster.is_empty():
            return ''
        else:
            if cluster.strand == '-':
                samflag = 16
            else:
                samflag = 0
            if not cluster.sequence:
                cluster.sequence = (len(cluster)+self.correction)*'A'
            if not cluster.is_singleton():
                lines = ''
                split_tags = cluster.absolute_split(0)
                for tag in split_tags:
                    tags = cluster.decompose()
                    for tag in tags:
                        lines = '%s%s'%(lines, sam_blueprint%(tag.name, samflag, tag.chromosome, tag.start+self.correction, len(tag), tag.sequence, (len(tag)+self.correction)*'?'))
                return lines
            else:
                return sam_blueprint%(cluster.name, samflag, cluster.chromosome, cluster.start+self.correction, len(cluster), cluster.sequence, (len(cluster)+self.correction)*'?')


class WigWriter(Writer):
    def write_line(self, cluster):
        lines = ''
        start = cluster.start
        for length, height in cluster:
            end = start+length-1
            if cluster.rounding:
                height = round(height)

            if height > 0:
                if cluster.rounding:
                    lines = '%s%s\t%s\t%s\t%.0f\n'%(lines,cluster.chromosome,start+self.correction, end, height)
                else:
                    lines = '%s%s\t%s\t%s\t%.2f\n'%(lines,cluster.chromosome,start+self.correction, end, height)
            start = end+1
        return lines

class FixedWigWriter(Writer):
    def write_line(self, cluster):
        lines = ''
        for length, height in cluster:
            if cluster.rounding:
                height = int(round(height))
            if height > 0:
                new_lines = ('%.2f\n'%height)*length
                lines = '%s%s'%(lines, new_lines)
        return lines

class VariableWigWriter(Writer):
    def write_line(self, cluster):
        """With this format we lose precision, but it loads faster in UCSC"""
        lines = ''
        height_array = []
        start = cluster.start
        for length, height in cluster:
            for i in range(length):
                height_array.append(height)

            while len(height_array) > self.span:
                mean_height = sum(height_array[0:self.span])/self.span
                for i in range(0, self.span):
                    height_array.pop(0)
                new_line = '%s\t%s\n'%(start, mean_height)
                lines = '%s%s'%(lines, new_line)
                start += self.span

        return lines

class PkWriter(Writer):
    def _format_line(self, cluster, start, acum_length, profile, max_height, max_height_pos, area):
        format_str = '%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s'%(cluster.chromosome, start+self.correction, start+acum_length-1, profile, max_height, cluster.strand, max_height_pos, area)
        if cluster.p_value != None: #it can be 0 and we want to print this, so DONT change this to "if cluster.p_value":
            format_str = '%s\t%s\n'%(format_str, cluster.p_value)
        else:
            format_str = '%s\n'%format_str
        return format_str


    def write_line(self, cluster):
        """prints a pk line out of the data in the object
        Also calculates the max height, max height position and area for all subclusters, replicating the functionality of max_height(), max_height_pos() and area()"""
        profile = ''
        lines = ''
        first = True
        acum_length = 0
        start = cluster.start
        max_height = -sys.maxint
        area = 0
        first_max_pos = 0
        last_max_pos = 0
        for length, height in cluster:
            if height > 0:   
                area += length*height           
                if first:
                    if cluster.rounding:
                        profile = '%s%s:%.0f'%(profile, length, height)
                    else:
                        profile = '%s%s:%.2f'%(profile, length, height)
                else:
                    if cluster.rounding:
                        profile = '%s|%s:%.0f'%(profile, length, height)
                    else:
                        profile = '%s|%s:%.2f'%(profile, length, height)

                acum_length += length
                if height > max_height:
                    max_height = height
                    first_max_pos = start + acum_length - length/2 - 1
                
                if height == max_height:
                    last_max_pos = start + acum_length - length/2 - 1

                first = False

            elif not first:
                lines = '%s%s'%(lines, self._format_line(cluster, start, acum_length, profile, max_height, (first_max_pos + last_max_pos)/2, area))
                start = start+acum_length+length
                acum_length = 0
                area = 0
                max_height = -sys.maxint
                first_max_pos = 0
                last_max_pos = 0
                first = True
                profile = ""

        if not first:
            lines = '%s%s'%(lines, self._format_line(cluster, start, acum_length, profile, max_height, (first_max_pos + last_max_pos)/2, area))

        return lines

######################################
#    CLUSTER OBJECT                  #
######################################


class Cluster(AbstractCore):
    """
    Represents one cluster of overlaping tags or a single Tag. This object can read and write in every format provided.
    It can be added, compared, subtracted to other cluster objects. Several interesting properties can be extracted from it.
    """
    def __init__(self, chromosome='', start=0, end=-1, strand='.', name='noname', score=0, rounding = False,
                 read=PK, write=PK, read_half_open=False, write_half_open=False, normalize_factor=1., tag_length=0, sequence=None, span=20, cached=False, verbose=False):
        """If you want the object to operate with integers instead
        of floats, set the rounding flag to True."""
        self.verbose = verbose
        self.chromosome = chromosome
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self._levels = []
        length = self.end-self.start+1
        if length > 0:
            self.add_level(length, 1)

        self.strand = strand
        self.rounding = rounding
        self.normalize_factor = normalize_factor
        self.score = score
        self.cached = cached
        self.read_as(read, read_half_open, cached)
        self.write_as(write, write_half_open, span)
        self.tag_length = tag_length
        self.sequence = sequence
        self._tag_cache = []
        self.read_count = 0
        self.strand = strand
        self.p_value = None

    def is_singleton(self):
        if self._tag_cache:
            self._flush_tag_cache()

        if self.is_empty():
            return False

        elif len(self._levels) > 1:
            return False

        elif self._levels[0][1] > 1:
            return False


        else:
            return True

    def read_as(self, format, half_open=False, cached=False):
        f = ReaderFactory()
        self.reader = f.create_reader(format, half_open, cached)

    def write_as(self, format, half_open=False, span=0):
        f = WriterFactory()
        self.writer = f.create_writer(format, half_open, span)

    def clear(self):
        """Empties the cluster"""
        self.chromosome = ''
        self.start = 0
        self.end = -1
        self.name = 'noname'
        del self._levels
        self._levels = []
        del self._tag_cache
        self._tag_cache = []
        self.strand = '.'
        self.score = 0
        self.read_count = 0

    def __len__(self):
        return self.end-self.start+1

    def __eq__(self, other):
        if  (self.chromosome != other.chromosome) or (self.start != other.start) or (self.strand != other.strand) or (self.end != other.end) or (self.name !=other.name) or (len(self._levels) != len(other._levels)):
            return False
        for i in xrange (0, len(self._levels)):
                if self._levels[i][0] != other._levels[i][0] or self._levels[i][1] != other._levels[i][1]:
                    return False
        return True

    def copy_cluster(self, copy_readers=True):
            """Returns a copy of the self cluster. Faster than copy.deepcopy()"""
            ret_cluster = Empty()
            ret_cluster.__class__ = self.__class__
            ret_cluster.start = self.start
            ret_cluster.end = self.end
            ret_cluster._levels = []
            ret_cluster._tag_cache = []

            if self._tag_cache: 
                for tag in self._tag_cache:
                    ret_cluster._tag_cache.append([tag[0], tag[1]])
            else:
                for length, height in self:
                    ret_cluster.add_level(length, height)

            ret_cluster.cached = self.cached
            ret_cluster.chromosome = self.chromosome
            if copy_readers:
                ret_cluster.read_as(self.reader.format, self.reader.half_open, self.cached)
                ret_cluster.write_as(self.writer.format, self.writer.half_open, self.writer.span)

            ret_cluster.score = self.score
            ret_cluster.name = self.name
            ret_cluster.strand = self.strand
            ret_cluster.rounding = self.rounding
            ret_cluster.normalize_factor = self.normalize_factor
            ret_cluster.tag_length = self.tag_length
            ret_cluster.sequence = self.sequence
            ret_cluster.read_count = self.read_count
            ret_cluster.verbose = self.verbose
            ret_cluster.p_value = self.p_value
            return ret_cluster

    #TODO raise 'there is no strand' error
    #TODO here I assume that the reads are shorter than the extension number.
    def extend(self, extension):
        if extension > 0:
            try:
                self._flush_tag_cache()
                if len(self._levels) > 1:
                    print 'Complex wig or pk clusters extension are not supported yet. Convert your files to bed format.'
                    raise InvalidLine

                if self.strand is '-':
                    previous_start = self.start
                    self.start = self.end - extension + 1
                    self._levels[0][0] += previous_start - self.start

                    if self.start < 1:
                        self.clear()
                        if self.verbose:
                            print 'Extending the line invalidates the read. Discarding ', self.chromosome, self.start, self.end
                        
                else:
                    previous_end = self.end
                    self.end = self.start + extension - 1
                    self._levels[-1][0] += self.end - previous_end

            except IndexError:
                #print 'Estoy atascado', self.start, self.end, self._levels
                pass

    def _subtrim(self, threshold, end, left):
        while len(self._levels) > 0:
            if self._levels[end][1] < threshold:
                if left:
                    self.start+= self._levels[end][0]
                else:
                    self.end -= self._levels[end][0]
                self._levels.pop(end)

            else:
                break
    
    def trim(self, threshold):
        """Trims the cluster to a given threshold"""
        self._subtrim(threshold, 0, True) #trim the left side of the cluster
        self._subtrim(threshold, -1, False) #trim the right side of the cluster

    def split(self, percentage=0.05, absolute=0):
        """
        Scans each cluster position from start to end and looks for local maxima x and local minima y.
        Given two consecutive local maxima x_{i} and x_{i+1} we define the smallest of them as x_{min}.
        For every y_{j} between two local maxima, the condition for splitting the cluster at y is defined as follows:

        y_{j}\leq x_{min}(1-t)

        Where t is a proportion threshold between 0 and 1. By default t=0.05. The cluster will divid at the local minimum. 
        """
        prev_height = -sys.maxint
        prev_local_maxima = 0
        new_cluster = self.copy_cluster()
        new_cluster.clear()
        new_cluster.start = self.start
        new_cluster.chromosome = self.chromosome
        clusters = []
        split_points = []
        going_up = True
        minimum_height = sys.maxint
        #get the local maxima
        level_number = 1
        for length, height in self:
            if height < prev_height: # we are going down
                if going_up: #if we were going up previously, we found a local maximum
                    if prev_local_maxima:
                        if absolute:
                            local_threshold = absolute
                        else:
                            local_threshold = min(prev_height, prev_local_maxima)*(1-percentage)
                        if minimum_height < local_threshold: #split point found
                            split_points.append(minimum_pos)               

                    prev_local_maxima = prev_height
                going_up = False
            else:
                if not going_up: #if we were going down previously, we reached a local minimum
                    minimum_height = prev_height
                    minimum_pos = level_number-1
                going_up = True

            prev_height = height
            prev_length = length
            level_number+=1
    
        if going_up: #if the cluster ended going up, it ended in a local maximum that was not detected previously
            if absolute:
                local_threshold = absolute
            else:
                local_threshold = min(prev_height, prev_local_maxima)*(1-percentage)
            if minimum_height < local_threshold: #split point found
                split_points.append(minimum_pos)    

        nucleotides = self.start
        if split_points:
            level_number = 1
            for length, height in self:
                if level_number in split_points:
                    right_len = length/2
                    left_len = right_len
                    if length%2 == 0: #its even
                        left_len-=1
                    if left_len:
                        new_cluster.add_level(left_len, height)
                    self.__subsplit(new_cluster, clusters, nucleotides, left_len+1)
                    if right_len:
                        new_cluster.add_level(right_len, height)
                else:
                    new_cluster.add_level(length, height) # add significant parts to the profile
                nucleotides+=length
                level_number+=1

            if not new_cluster.is_empty():
               self.__subsplit(new_cluster, clusters, nucleotides, length)
            return clusters
        else:
            return [self]



        return clusters

    def __subsplit(self, new_cluster, clusters, nucleotides, length):
        """sub method of split"""
        new_cluster.chromosome=self.chromosome
        new_cluster._recalculate_end()
        clusters.append(new_cluster.copy_cluster())
        new_cluster.clear()
        new_cluster.start=nucleotides+length

    def absolute_split(self, threshold):
        """Returns the original cluster or several clusters if we find subclusters"""
        if threshold is None:
            threshold=float(self.max_height())*percentage
            if threshold < 1: #or at least 1 nucleotide
                threshold = 1

        new_cluster = self.copy_cluster()
        new_cluster.clear()
        new_cluster.start = self.start
        nucleotides = self.start
        clusters = []
        for length, height in self:
            if height<=threshold:
                if new_cluster.is_empty():
                    new_cluster.start=nucleotides+length # trim insignificant beginning of a cluster
                else:
                    self.__subsplit(new_cluster, clusters, nucleotides, length)
            else:
                new_cluster.add_level(length, height) # add significant parts to the profile
            nucleotides+=length

        if not new_cluster.is_empty():
           self.__subsplit(new_cluster, clusters, nucleotides, length)

        return clusters

    def is_artifact(self, condition = 0.3):
        """Returns True if the cluster is considered an artifact.
           It is considered an artifact if its shorter than 100 nucleotides,
           or the maximum height length is more than 30% of the cluster (block cluster)"""
        self._flush_tag_cache()
        if len(self) < 100:
            return True
        h = self.max_height()
        max_len = 0.
        for length, height in self:
            if h == height:
                max_len = max(float(length), max_len)

        return max_len/len(self) > condition

    def is_significant(self, threshold, poisson_type="height"):
        """Returns True if the cluster is significant provided a threshold, otherwise False"""
        if poisson_type=="height":
            return threshold <= int(round(self.max_height()))
        else:
            return threshold <= int(round(self.area()))

    def read_line(self, line):
        self.reader.read_line(self, line)

    def write_line(self):
        if self._tag_cache:
            self._flush_tag_cache()
        return self.writer.write_line(self)

    def set_tag_length(self, tag_length):
        self.tag_length = tag_length

    def decompose(self):
        """Returns a list of Singletons that form the cluster given a tag length. If the tag length is not provided the cluster was formed by tags with
        different lengths, its impossible to reconstruct the tags due to ambiguity (multiple different combination of tags can be used to obtain the same cluster)"""
        if self.tag_length is not None:
            tags = []
            cluster_copy = self.copy_cluster()
            dummy_cluster = self.copy_cluster()
            dummy_cluster._levels = []
            dummy_cluster.add_level(self.tag_length, 1)
            dummy_cluster.start = cluster_copy.start
            dummy_cluster.end = cluster_copy.start+self.tag_length-1
            while cluster_copy.max_height() > 0:
                cluster_copy -= dummy_cluster
                tags.append(dummy_cluster.copy_cluster())
                dummy_cluster.start = cluster_copy.start
                dummy_cluster.end = cluster_copy.start+self.tag_length-1

        else:
            raise InsufficientData

        return tags

    def __sub__(self, other):
        if debug: print "Flushing experiment..."
        if self._tag_cache:
            self._flush_tag_cache()
        if debug: print "Flushing control..."
        if other._tag_cache:
            other._flush_tag_cache()
        if debug: print "Copying..."
        result = self.copy_cluster()
        if debug: print "Subtracting..."
        if result.chromosome == other.chromosome:
            other_acum_levels = 0
            j = 0
            while (len(other._levels) > j):
                if debug: print "j:", j, len(other._levels)
                other_level_start = other.start + other_acum_levels
                other_acum_levels += other._levels[j][0] #for speed
                other_level_end = other.start + other_acum_levels
                other_height = other._levels[j][1] #for speed
                other_length = other._levels[j][0] #for speed
                i = 0
                acum_levels = 0
                while (len(result._levels) > i):
                    level_start = result.start + acum_levels
                    acum_levels += result._levels[i][0]
                    level_end = result.start + acum_levels
                    height = result._levels[i][1]

                    if (other_level_start >= level_start and other_level_start < level_end) or (other_level_end > level_start and other_level_end < level_end) or (other_level_start <= level_start and other_level_end >= level_end):

                        if other_level_start <= level_start and other_level_end >= level_end:
                            result._levels[i][1] -= other_height

                        elif other_level_start <= level_start and other_level_end < level_end:
                            result._levels[i][0] = other_level_end - level_start
                            result._levels[i][1] -= other_height
                            result._levels.insert(i+1, [level_end - other_level_end, height])
                            level_end = other_level_end

                        elif other_level_start > level_start and other_level_start < level_end and other_level_end >= level_end:
                            result._levels[i][0] -= level_end-other_level_start
                            result._levels.insert(i+1, [level_end-other_level_start,result._levels[i][1]-other_height])
                            i+=1

                        elif other_level_start > level_start and other_level_end < level_end:
                            result._levels[i][0] = other_length
                            result._levels[i][1] -= other_height
                            result._levels.insert(i+1, [level_end - other_level_end, height])
                            result._levels.insert(i, [other_level_start - level_start, height])

                    i+=1

                j+=1
            if debug: print "Cleaning..."
            result._clean_levels()
        return result

    def __iter__(self):
        if self._tag_cache:
            self._flush_tag_cache()
        return self._levels.__iter__()


    def max_height(self):
        """Returns the maximum height in the cluster"""
        max_height = 0.
        for length, height in self:
            max_height = max(max_height, height)
        return max_height


    def max_height_pos(self):
        """
        Returns the position where the maximum height is located.
        The central positions of the first and the last maximum are calculated.
        The max_height_pos will be in the middle of these two positions.
        """
        max_height = 0
        acum_length = 0
        first_max = 0
        last_max = 0
        for length, height in self:
            acum_length += length
            if round(height) > round(max_height):
                max_height = height
                first_max = self.start + acum_length - length/2 - 1

            if round(height) == round(max_height):	
                last_max = self.start + acum_length - length/2 - 1

        pos = (first_max + last_max)/2
        return pos 
 

    def area(self):
        """
        Returns the area of the peak
        """

        sum_area = 0
        for length, height in self:
            sum_area += length*height

        return sum_area

    def add_level(self, length, height):
        self._levels.append([int(length), float(height)])


    def __add__(self, other):
        """
        Adds the levels of 2 selfs, activating the + operator for this type of object results = self + self2
        """
        if other._tag_cache:
            self._tag_cache.extend(other._tag_cache)
            return self
        else: 
            result = self.copy_cluster()
            if result.read_count and other.read_count:
                result.read_count += other.read_count #if both clusters have a read count, add it
            other_acum_levels = 0
            #add zeros so both selfs have equal length and every level is added
            if other.end > result.end:
                result.add_level(other.end-result.end, 0)

            if other.start < result.start:
                result._levels.insert(0, [result.start - other.start, 0])
                result.start = other.start

            for j in xrange(0, len(other._levels)):
                other_level_start = other.start + other_acum_levels
                other_acum_levels += other._levels[j][0]
                other_level_end = other.start + other_acum_levels
                other_height = other._levels[j][1]
                other_length = other._levels[j][0]
                i = 0
                acum_levels = 0
                while (len(result._levels) > i):
                    level_start = result.start + acum_levels
                    acum_levels += result._levels[i][0]
                    level_end = result.start + acum_levels
                    height = result._levels[i][1]
                    if (other_level_start >= level_start and other_level_start < level_end)  or (other_level_end > level_start and other_level_end < level_end) or (other_level_start <= level_start and other_level_end >= level_end):
                        if other_level_start <= level_start and other_level_end >= level_end:
                            result._levels[i][1] += other_height

                        elif other_level_start <= level_start and other_level_end < level_end:
                            result._levels[i][0] = other_level_end - level_start
                            result._levels[i][1] += other_height
                            result._levels.insert(i+1, [level_end - other_level_end, height])
                            level_end = other_level_end

                        elif other_level_start > level_start and other_level_start < level_end and other_level_end >= level_end:
                            result._levels[i][0] -= level_end-other_level_start
                            result._levels.insert(i+1, [level_end-other_level_start,result._levels[i][1]+other_height])
                            i+=1

                        elif other_level_start > level_start and other_level_end < level_end:
                            result._levels[i][0] = other_length
                            result._levels[i][1] += other_height
                            result._levels.insert(i+1, [level_end - other_level_end, height])
                            result._levels.insert(i, [other_level_start - level_start, height])
                    i+=1

            result._clean_levels()
        return result

    def _flush_tag_cache(self):
        """
        Joins all reads in levels. Assumes that the tags were sorted (TODO: Is there a FAST way of detecting if they are sorted, so we could sort them if neccesary?)
        """
        import heapq
        array_ends = []
        previous_start = -1
        smallest_end = sys.maxint
        for current_start, current_end in self._tag_cache:
            while current_start > smallest_end and len(array_ends) > 0:
                self._levels.append([smallest_end-previous_start+1, len(array_ends)*self.normalize_factor])
                previous_start = heapq.heappop(array_ends)
                try:
                    while previous_start == heapq.nsmallest(1, array_ends)[0]:
                        previous_start = heapq.heappop(array_ends)
                except IndexError:
                    pass #if array_ends is empty, go on
                previous_start = previous_start + 1
                if len(array_ends) > 0:
                    smallest_end = heapq.nsmallest(1, array_ends)[0]

            if previous_start != current_start and len(array_ends) > 0:
                self._levels.append([current_start-previous_start, len(array_ends)*self.normalize_factor])
            previous_start = current_start
            heapq.heappush(array_ends, current_end)
            smallest_end = heapq.nsmallest(1, array_ends)[0]

        if len(array_ends) > 0:
            previous_end = -1
            while len(array_ends) > 0:
                if previous_end != smallest_end and len(array_ends) > 0:
                    self._levels.append([smallest_end-previous_start+1, len(array_ends)*self.normalize_factor])
                previous_start = heapq.heappop(array_ends)+1
                if len(array_ends) > 0:
                    previous_end = smallest_end
                    smallest_end = heapq.nsmallest(1, array_ends)[0]

        self._tag_cache = []
        self._clean_levels()
        

    def get_heights(self):
        """returns all the heights in an array"""
        ret = []
        for length, height in self:
            for i in xrange(0, length):
                ret.append(height)

        return ret

    def is_empty(self):
        """Returns True if the Cluster object is empty, returns False otherwise."""
        return len(self._levels) == 0 and len(self._tag_cache) == 0

    def _clean_levels(self):
        """
        Cleans the levels array joining the values that have the same height and discarding the values
        that have height < 0 at the initial positions.
        """
        previous_height = -1
        previous_length = 0
        i = 0
        if len(self._levels) > 0:
            while self._levels[0][0] <= 0:
                self._levels.pop(0)
                if self.is_empty():
                    break



        if len(self._levels) > 0:
            while self._levels[0][1] <= 0: #delete the 0 to the left of the Cluster
                self.start += self._levels[0][0]
                self._levels.pop(0)
                if self.is_empty():
                    break
        

        if len(self._levels) > 0:
            while self._levels[-1][1] <= 0: #delete the 0 to the right of the Cluster
                self._levels.pop(len(self._levels)-1)
                if len(self._levels) == 0:
                    break

        while (len(self._levels) > i): #join all identical levels
            if self._levels[i][1] == previous_height:
                self._levels[i][0] += previous_length
                previous_length = self._levels[i][0]
                self._levels.pop(i-1)

            else:
                previous_length = self._levels[i][0]
                previous_height = self._levels[i][1]
                i+=1


        self._recalculate_end()

    def _recalculate_end(self):
        """Recalculates the end value to avoid possible incosistency in the data"""
        self.end = self.start-1
        for length, height in self:
            self.end += length

    def has_duplicates(self, limit=20):
        """Returns true if the cluster has any ocurrence that exceeds the number of duplicated reads
        specified in the limit variable"""
        if limit < 1:
            return False
        else:
            previous_height = 0
            for length, height in self:
                if height > previous_height + limit:
                    return True
                previous_height = height
            return False

    def __str__(self):
        return "<Cluster object> chr: %s start: %s end: %s name: %s "%(self.chromosome, self.start, self.end, self.name)


#######################
#   REGION  OBJECT    #
#######################
class Region(AbstractCore):
    def __init__(self, start=0, end=0, chromosome=None, strand=None):
        self.start = int(start)
        self.end = int(end)
        self.chromosome = chromosome
        self.strand = strand
        self.tags = []
        self.clusters = []

    def rpkm(self, total_reads):
        """Original definition: Reads per kilobase of exon model per million mapped reads. We generalize to: Reads per kilobase of region per million mapped reads. Added 1 pseudocount per region to avoid 0s"""
        return (10e9*float(len(self)+1))/((len(self.tags)+1)*total_reads)


    def __sub_swap(self, region, swap1, swap2):
        for tag in region.tags:
            if random.randint(0,1):
                swap1.add_tags(tag)
            else:
                swap2.add_tags(tag)

    def swap(self, region_b):
        "Given 2 regions, return 2 new regions with the reads of both regions mixed aleatoriely"
        swap1 = Region(self.start, self.end, self.chromosome)
        swap2 = Region(self.start, self.end, self.chromosome)
        self.__sub_swap(self, swap1, swap2)
        self.__sub_swap(region_b, swap1, swap2)
        return (swap1, swap2)        

    def __str__(self):
        return "chr: %s start: %s end: %s number_of_tags: %s"%(self.chromosome, self.start, self.end, len(self.tags))

    def _numpos_higher_than(self, h, nis):
        ret = 0
        for key, value in nis.items():
            if key >= h:
              ret += value
        return float(ret)

    def _get_probability(self, nis, h):
        return self._numpos_higher_than(h, nis)/self._numpos_higher_than(1, nis)

    def __len__(self):
        return self.end-self.start+1

    def copy(self):
        """Copies a Region object into another"""
        new_region = Empty()
        new_region.__class__ = self.__class__
        new_region.start = self.start
        new_region.end = self.end
        new_region.chromosome = self.chromosome
        new_region.tags = []
        new_region.clusters = []
        new_region.add_tags(self.tags)
        new_region.strand = self.strand
        for cluster in self.clusters:
            new_region.clusters.append(cluster)
        return new_region

    def num_tags(self):
        return len(self.tags)

    def max_height(self):
        max_height = 0
        for cluster in self.clusters:
            max_height = max(max_height, cluster.max_height())
        return max_height

    def tag_len_average(self):
        ret = 0.
        for tag in self.tags:
            ret += len(tag)
        return ret/len(self.tags)            

    def binomial(self, n, p, k): 
        return (math.factorial(n)*p**k*(1-p)**(n-k))/(math.factorial(k)*math.factorial(n-k))

         
    def p0(self, l, N):
        return (N-l)/N

    def p1(self, l, N):
        return l/N
    
    def join(self, other):
        """Joins two regions. Works with a cluster object too (if flushed properly)"""
        self.start = min(self.start, other.start)
        self.end = max(self.end, other.end)
   
    def get_FDR_clusters(self, repeats=100, masker_tags=[]):

        max_height = int(self.max_height())
        #Get the N values and the probabilities for the real and for the simulated random instances
        real_heights = self.get_heights()
        #we do the same calculation but shuffling the tags randomly. We do it as many times as the variable repeats asks for
        Pr_of_h = defaultdict(int)
        random_variances = defaultdict(int)
        random_region = Region(self.start, self.end)
        random_region.add_tags(self.tags)
        #Get the repeat regions that overlap with the region
        """
        masker_tags = repeat_reader.get_overlaping_clusters(region, overlap=0.000001) #get all overlaping masker tags
        masker_region = Region(region.start, region.end, region.chromosome)
        masker_region.add_tags(masker_tags)
        """

        for r in xrange(0, repeats):
            sys.stdout.flush()
            random_region.shuffle_tags(masker_tags)
            nis = random_region.get_heights()
            for h in xrange(1, max_height+1):
                if not random_variances[h]:
                    random_variances[h] = []
                prob = self._get_probability(nis, h)
                Pr_of_h[h] += prob
                random_variances[h].append(prob)


        #Get the mean and the variance for the random iterations and calculate the FDR
        found = False
        significant_height = sys.maxint
        p_values_per_height = defaultdict(int) #if it doesnt have the values, will return a 0

        """
        number_of_reads = len(self.tags)
        average_tag_length = self.tag_len_average()
        region_length = len(self)
         
        est_pr_of_h = defaultdict(int)
        acum = 1
        for h in xrange(1, max_height+1):
            est_pr_of_h[h] = acum
            acum -= self.binomial(number_of_reads, (average_tag_length/region_length), h)
        """

        for h in xrange(1, max_height+1):
            random_mean = Pr_of_h[h]/repeats
            #print "random_prob:", h, random_mean
            #print "Binomial:", h, est_pr_of_h[h]
            #print 

            random_variance = 0
            h_prob = self._get_probability(real_heights, h)
            for i in xrange(0, repeats):
                random_variance += (random_variances[h][i]-random_mean)**2
            random_variance = math.sqrt(random_variance/repeats)
            p_value = random_mean+random_variance/h_prob
            if p_value < 0.00000001: #A smaller p-value than 10e-8 is irrelevant. For speed.
                break
            p_values_per_height[h] = random_mean+random_variance/h_prob
        #add the p_values to the clusters and return them
        ret_clusters = self.get_clusters()
        for i in range(0 , len(ret_clusters)):
            ret_clusters[i].p_value = p_values_per_height[int(ret_clusters[i].max_height())]

        return ret_clusters

    def get_heights(self):
        """Get the number of nucleotides that have a certain height. Returns a defaultdict"""
        heights_dict = defaultdict(int)
        for cluster in self.clusters:
            cluster._flush_tag_cache()
            for length, height in cluster:
                heights_dict[int(height)] += length
        return heights_dict


    def shuffle_tags(self, masker_tags=[]):
        """Shuffle the present tags in the region randomly"""
        for tag in self.tags:
            repeat = True
            repeat_count = 0
            #length = self.end-self.start #self._recalculate_end() is too slow
            while (repeat):
                #Try to find a random spot that doesnt fall inside a repeat region.
                # If a repeat region is hit more than 2000 times, surrender and leave the tag where it last was randomized
                #(almost impossible with normal genomes, but it may be possible if a gene with more than
                #99.9% of repeats is analyzed) 0.999**2000=0.13519992539749945  0.99**2000=1.8637566029922332e-09
                repeat = False
                tag.start = random.randint(self.start, self.end)
                tag._recalculate_end()
                for repeat_tag in masker_tags:
                    if tag.overlap(repeat_tag) > 0:
                        repeat = True
                        repeat_count+=1
                        break

                if repeat_count > 2000:
                    print 'Couldnt find a suitable randomization spot after 2000 tries, surrendering'
                    break

        #recreate the clusters
        self.clusterize()

    def _sub_add_tag(self, tag):
        if not self.strand or tag.strand == self.strand:
            self.tags.append(tag)
        

    def add_tags(self, tags, clusterize=False):
        """This method reads a list of tags or a single tag (Cluster objects, not unprocessed lines). If strand is set, then only the tags with the selected strand are added"""
        if type(tags)==type(list()):
            for tag in tags:
                self._sub_add_tag(tag)
        elif type(tags)==type(Cluster()):
            self._sub_add_tag(tags)
        else:
            print 'Invalid tag. Tags need to be either Cluster or List objects'

        if clusterize:
            self.clusterize()

    def clusterize(self):
        """Creates the Cluster objects of the tags in the Region object"""
        self.clusters = []
        self.tags.sort(key=lambda x: (x.start, x.end))
        if self.tags:
            #Insert first cluster object
            self.clusters.append(Cluster(read=self.tags[0].reader.format, cached=True))
            for i in xrange(0, len(self.tags)):
                if not self.tags[i].is_empty():
                    if not self.clusters[-1].overlap(self.tags[i]) > 0 and not self.clusters[-1].is_empty():
                        self.clusters.append(Cluster(read=self.tags[0].reader.format, cached=True)) #create a new cluster object
                    try:
                        self.clusters[-1].read_line(self.tags[i].write_line())
                    except InvalidLine:
                        print "A VER:", self.tags[i].write_line()
                        raise


    def percentage_covered(self):
        """Returns the percentage of the region covered by tags"""
        covered = 0.
        total = len(self)
        for cluster in self.clusters:
            for length, height in cluster:
                covered+=length
        if total > 0:
            return covered/total
        else:
            return 0


    def get_metacluster(self):
        """Returns a cluster object that contains the levels of all previous clusters combined, with gaps (zeros) between them"""
        if self.clusters:
            if self.clusters[0]._tag_cache:
                self.clusters[0]._flush_tag_cache()
            metacluster = self.clusters[0].copy_cluster()
            previous_end = metacluster.end
            for i in range(1, len(self.clusters)):
                self.clusters[i]._flush_tag_cache() #Need the end in its proper place
                if self.clusters[i].start > previous_end: # add zeros between levels
                    metacluster.add_level(self.clusters[i].start-previous_end-1, 0)
            
                for length, height in self.clusters[i]:
                    metacluster.add_level(length, height)
                previous_end = self.clusters[i].end

            metacluster._recalculate_end()
            return metacluster
        else:
            return Cluster()
        

    def get_clusters(self, height=1):
        """Gets the clusters inside the region higher than the marked height. By default returns all clusters with at least 2 reads"""
        significant_clusters = []
        for cluster in self.clusters:
            if cluster.max_height() > height:
                significant_clusters.append(cluster)
        return significant_clusters

    def write(self):
        """Returns a line in bed format"""
        strand = self.strand
        if not strand:
            strand = '.'
        return "%s\t%s\t%s\t0\tpyicos_region\t%s\n"%(self.chromosome, self.start, self.end, strand)


