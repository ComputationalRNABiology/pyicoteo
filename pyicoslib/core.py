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


debug = False
verbose = False

#MIRO = 'miro'
ELAND = 'eland'
ELAND_EXPORT = 'eland_export'
BED = 'bed'
WIG = 'bed_wig'
VARIABLE_WIG = 'variable_wig'
FIXED_WIG = 'fixed_wig'
PK = 'bedpk'
SPK = 'bedspk'

CLUSTER_FORMATS = (WIG, VARIABLE_WIG, FIXED_WIG, PK, SPK)
WIG_FORMATS = (WIG, VARIABLE_WIG, FIXED_WIG)

READ_FORMATS = (ELAND, BED, WIG, PK, SPK) #formats that we actually can read as
WRITE_FORMATS = (ELAND, BED, WIG, VARIABLE_WIG, PK, SPK) #formats we can actually write as


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
            cluster.strand == ''

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
            cluster.strand = ''

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
                    self._add_strand(cluster, line)
                else:
                    cluster._tag_cache.append([new_start, int(line[2])])
                    cluster.start = min(cluster.start, new_start)
                    cluster.end = max(cluster.end, int(line[2]))
                    if len(line) > 5:
                        if cluster.strand != line[5]:
                            cluster.strand == ''
                    else:
                        cluster.strand == ''
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
                    print cluster.normalize_factor
                    cluster.add_level(cluster.end-cluster.start+1, cluster.normalize_factor)

        except ValueError:
            raise InvalidLine

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
        except ValueError:
            raise InvalidLine

class ElandReader(Reader):
    qualityFilter = re.compile(r'chr\w{1,2}.fa')
    
    #sort_func = lambda x:(x.split()[6], int(x.split()[7]), len(x.split()[1]))
    def read_line(self, cluster, line):
        if line is not None and line != '\n':
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
                    if line[8] is 'F':
                        cluster.strand = '+'
                    else:
                        cluster.strand = '-'

                else:
                    self._add_line_to_cluster(line, cluster)

            except ValueError:
                raise InvalidLine
            except IndexError:
                raise InvalidLine

    def eland_quality_filter(self, line):
        """checks if the eland line passes the quality conditions"""
        return self.qualityFilter.search(line)

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

                    if len(line) > 4:
                        cluster.strand = line[4]
                    else:
                        cluster.strand = ""

                    cluster._recalculate_end()
        except ValueError:
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
    def _format_line(self, cluster, start, acum_length, profile):
        if self.format == PK:
            return '%s\t%s\t%s\t%s\t%s\t.\t%s\t%s\n'%(cluster.chromosome, start+self.correction, start+acum_length-1, profile, cluster.get_max_height(), cluster.get_max_height_pos(), cluster.get_area())
        else: #Its SPK
            return '%s\t%s\t%s\t%s\t%s\t%s\n'%(cluster.chromosome, start+self.correction, start+acum_length-1, profile, cluster.score, cluster.strand)

    def write_line(self, cluster):
        """prints a pk line out of the data in the object"""
        profile = ''
        lines = ''
        first = True
        acum_length = 0
        start = cluster.start
        for length, height in cluster:
            if cluster.rounding:
                height = round(height)
            if height > 0:
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
                first = False

            elif not first:
                lines = '%s%s'%(lines, self._format_line(cluster, start, acum_length, profile))
                start = start+acum_length+length
                acum_length = 0
                first = True
                profile = ""

        if not first:
            lines = '%s%s'%(lines, self._format_line(cluster, start, acum_length, profile))

        return lines

######################################
#    CLUSTER OBJECT                  #
######################################


class Cluster:
    """
    Represents one cluster of overlaping tags or a single Tag. This object can read and write in every format provided.
    It can be added, compared, subtracted to other cluster objects. Several interesting properties can be extracted from it.
    """
    def __init__(self, chromosome='', start=0, end=-1, strand='', name='noname', score=0, rounding = False,
                 read=PK, write=PK, read_half_open=False, write_half_open=False, normalize_factor=1., tag_length=None, sequence=None, span=20, cached=False, verbose=False):
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
        self.read_as(read, read_half_open)
        self.write_as(write, write_half_open, span)

        self.tag_length = tag_length
        self.sequence = sequence
        self._tag_cache = []
        self.read_count = 0

    def is_singleton(self):
        if self._tag_cache:
            self._flush_tag_cache()
        if self.is_empty():
            return False

        elif len(self._levels) > 1:
            return False

        elif self._levels[0][1] > 1:
            return False

        elif self.tag_length:
            if self._levels[0][0] > self.tag_length-1:
                return False

        else:
            return True

    def read_as(self, format, half_open=False, cached=True):
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
        self._levels = []
        self._tag_cache = []
        self.strand = ''
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

    def copy_cluster(self):
            """Returns a copy of the self cluster. Faster than copy.deepcopy()"""
            ret_cluster = Empty()
            ret_cluster.__class__ = self.__class__
            ret_cluster.start = self.start
            ret_cluster.end = self.end
            ret_cluster._levels = []
            ret_cluster._tag_cache = []
            for length, height in self:
                ret_cluster.add_level(length, height)
            ret_cluster.chromosome = self.chromosome
            ret_cluster.read_as(self.reader.format, self.reader.half_open)
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

    def split(self, percentage=0.05):
        """
        Scans each cluster position from start to end and looks for local maxima x and local minima y.
        Given two consecutive local maxima x_{i} and x_{i+1} we define the smallest of them as x_{min}.
        For every y_{j} between two local maxima, the condition for splitting the cluster at y is defined as follows:

        y_{j}\leq x_{min}(1-t)

        Where t is a proportion threshold between 0 and 1. By default t=0.05, meaning that by default clusters will be split if the "gap"
        in between has a signal 5% smaller than the smaller local maxima.
        """
        prev_height = -sys.maxint
        prev_local_maxima = None
        new_cluster = self.copy_cluster()
        new_cluster.clear()
        new_cluster.start = self.start
        new_cluster.chromosome = self.chromosome
        clusters = []
        thresholds = []
        thresholds.append(-sys.maxint)
        local_threshold = sys.maxint
        going_up = True
        #get the local maxima
        for length, height in self:
            if height < prev_height: # we are going down
                if going_up: #if we were going up previously, we found a local maxima
                    if prev_local_maxima:
                        local_threshold = min(prev_height, prev_local_maxima)*(1-percentage)
                        thresholds.append(local_threshold)
                    prev_local_maxima = prev_height
                going_up = False
            else:
                going_up = True

            prev_height = height

        thresholds.append(-sys.maxint)
        prev_height = -sys.maxint
        maxima_count = 0
        nucleotides = self.start
        going_up = True
        #split using the local maxima threshold information
        for length, height in self:
            if height < prev_height: 
                if going_up: #previous was a maxima
                    maxima_count += 1
                going_up = False
            else:
                going_up = True

            if height <= thresholds[maxima_count]:
                if not new_cluster.is_empty():
                    self.__subsplit(new_cluster, clusters, nucleotides, length)
                new_cluster.start=nucleotides+length
            else:
                 new_cluster.add_level(length, height) # add significant parts to the profile

            prev_height = height
            nucleotides += length

        if not new_cluster.is_empty():
           self.__subsplit(new_cluster, clusters, nucleotides, length)

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
            threshold=float(self.get_max_height())*percentage
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
        if len(self) < 100:
            return True
        h = self.get_max_height()
        max_len = 0.
        for length, height in self:
            if h == height:
                max_len = max(float(length), max_len)

        return max_len/len(self) > condition

    def is_significant(self, threshold):
        """Returns True if the cluster is significant provided a threshold, otherwise False"""
        return threshold <= self.get_max_height()

    def intersects(self, other):
        """Returns true if a Cluster intersects with another Cluster"""
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
            while cluster_copy.get_max_height() > 0:
                cluster_copy -= dummy_cluster
                tags.append(dummy_cluster.copy_cluster())
                dummy_cluster.start = cluster_copy.start
                dummy_cluster.end = cluster_copy.start+self.tag_length-1

        else:
            raise InsufficientData

        return tags

    def __sub__(self, other):
        if other._tag_cache:
            other._flush_tag_cache()

        result = self.copy_cluster()
        if result.chromosome == other.chromosome:
            other_acum_levels = 0
            j = 0
            while (len(other._levels) > j):
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

            result._clean_levels()

        return result

    def __iter__(self):
        if self._tag_cache:
            self._flush_tag_cache()
        return self._levels.__iter__()

    def get_max_height(self):
        """Returns the maximum height in the cluster"""
        max_height = 0.
        for length, height in self:
            max_height = max(max_height, height)
        return max_height

    def get_max_height_pos(self):
	# changed by Sonja (to Juanra's idea)
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
            if int(height) > int(max_height):
                max_height = height
                first_max = self.start + acum_length - length/2 - 1
	    if int(height) == int(max_height):	
                last_max = self.start + acum_length - length/2 - 1

        pos = (first_max + last_max)/2
        return pos 
 

    def get_area(self):
	# added by Sonja
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
            other._flush_tag_cache()
            
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
        import heapq
        array_ends = []
        previous_start = -1
        smallest_end = sys.maxint
        for current_start, current_end in self._tag_cache:
            while current_start > smallest_end and len(array_ends) > 0:
                self._levels.append([smallest_end-previous_start+1, len(array_ends)*self.normalize_factor])
                previous_start = heapq.heappop(array_ends)
                while previous_start == heapq.nsmallest(1, array_ends)[0]:
                    previous_start = heapq.heappop(array_ends)
                
                previous_start = previous_start + 1
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


#######################
#   REGION  OBJECT    #
#######################
class Region:
    def __init__(self, start=0, end=0, name=None, chromosome=None):
        self.start = int(start)
        self.end = int(end)
        self.name = name
        self.chromosome = chromosome
        self.tags = []
        self.clusters = []

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

    def num_tags(self):
        return len(self.tags)

    def get_max_height(self):
        max_height = 0
        for cluster in self.clusters:
            max_height = max(max_height, cluster.get_max_height())
        return max_height

    def get_FDR_clusters(self, repeats=100, p_value=0.01, masker_tags=[]):
        significant_clusters = []
        max_height = int(self.get_max_height())
        #Get the N values and the probabilities for the real and for the simulated random instances
        real_heights = self.get_heights()
        #we do the same calculation but shuffling the tags randomly. We do it as many times as the variable repeats asks for
        Pr_of_h = defaultdict(int)
        random_variances = defaultdict(int)
        random_region = Region(self.start, self.end)
        random_region.add_tags(self.tags)
        #Get the repeat regions that overlap with the region
        """masker_tags = repeat_reader.get_overlaping_clusters(region, overlap=0.000001) #get all overlaping masker tags
        masker_region = Region(region.chromosome, region.start, region.end, region.name)
        masker_region.add_tags(masker_tags)"""
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
        for h in xrange(1, max_height+1):
            random_mean = Pr_of_h[h]/repeats
            random_variance = 0
            h_prob = self._get_probability(real_heights, h)
            for i in xrange(0, repeats):
                random_variance += (random_variances[h][i]-random_mean)**2
            random_variance = math.sqrt(random_variance/repeats)
            FDR = random_mean+random_variance/h_prob
            if FDR < p_value:
                if not found: #If a p-value smaller than the target is reached, stop the calculations and print the result
                    significant_height = h
                    found = True
                    break    
        for cluster in self.get_clusters():
            if cluster.get_max_height() > significant_height:
                significant_clusters.append(cluster)
        return significant_clusters

    def get_heights(self):
        """Get the number of nucleotides that have a certain height. Returns a defaultdict"""
        heights_dict = defaultdict(int)
        for cluster in self.clusters:
            cluster._flush_tag_cache()
            for length, height in cluster:
                heights_dict[int(height)] += length
        return heights_dict

    def shuffle_tags(self, masker_tags=[]):
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
                    print 'I Surrender!'
                    break

        #we need to sort the tags
        self._create_clusters()

    def add_tags(self, tags):
        """This method reads a list of tags or a single tag"""
        if type(tags)==type(list()):
            self.tags.extend(tags)
        elif type(tags)==type(Cluster()):
            self.tags.append(tags)
        else:
            print 'Invalid tag!!!'
        self._create_clusters()

    def _create_clusters(self):
        """Creates the Cluster objects that will hold the profile of the region."""
        self.clusters = []
        self.tags.sort(key=lambda x: (x.start, x.end))
        if self.tags:
            #get the first tag
            self.clusters.append(Cluster(read=self.tags[0].reader.format, cached=True))

            for i in xrange(0, len(self.tags)):
                if self.clusters[-1].overlap(self.tags[i]) > 0 or self.clusters[-1].is_empty():
                    self.clusters[-1].read_line(self.tags[i].write_line())
                else:
                    self.clusters.append(Cluster(read=self.tags[0].reader.format, cached=True))
    """
    def _create_clusters(self):
        self.clusters = []
        self.tags.sort(key=lambda x: (x.start, x.end))
        for tag in self.tags:
            self.clusters.append(tag)
        i = 0
        while i < len(self.clusters)-1:
            if self.clusters[i].intersects(self.clusters[i+1]) or self.clusters[i].intersects(self.clusters[i+1]) or self.clusters[i+1].intersects(self.clusters[i]):
                self.clusters[i] += self.clusters[i+1]
                self.clusters.pop(i+1)
            else:
                i += 1
    """


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


    def get_complete_profile(self):
        """Incomplete: Needs 0s between clusters"""
        profile = ''
        previous_end = -1
        for cluster in self.clusters:
            for length, height in cluster:
                profile = '%s|%s:%s'%(profile, length, height)
        return profile[1:].strip()

    def get_clusters(self, height=1):
        """Gets the clusters inside the region higher than the height"""
        significant_clusters = []
        for cluster in self.clusters:
            if cluster.get_max_height() >= height:
                significant_clusters.append(cluster)
        return significant_clusters
