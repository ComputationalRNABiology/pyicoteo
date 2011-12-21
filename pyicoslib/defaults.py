#Formats
ELAND = 'eland'
ELAND_EXPORT = 'eland_export'
BED = 'bed'
BED12 = 'bed12'
WIG = 'bed_wig'
VARIABLE_WIG = 'variable_wig'
FIXED_WIG = 'fixed_wig'
PK = 'bed_pk'
SPK = 'bed_spk'
SAM = 'sam'
COUNTS = 'counts'
CLUSTER_FORMATS = (WIG, VARIABLE_WIG, FIXED_WIG, PK, SPK)
WIG_FORMATS = (WIG, VARIABLE_WIG, FIXED_WIG)
READ_FORMATS = (ELAND, BED, WIG, PK, SPK, SAM, COUNTS) #formats that we actually can read as
WRITE_FORMATS = (ELAND, BED, WIG, VARIABLE_WIG, PK, SPK) #formats we can actually write as

REGION_FORMATS = (BED, BED12)
#Enrichment header
enrichment_keys  = ['name', 'start', 'end', 'name2', 'score', 'strand', 'signal_a', 'signal_b', 'signal_prime_1', 'signal_prime_2',
                    'A','M','total_reads_a','total_reads_b','num_tags_a','num_tags_b','A_prime','M_prime',
                    'total_reads_a','total_reads_b','total_reads_background_1','total_reads_background_2', 'A_median', 'mean', 'sd', 'zscore']

#Default values for parser
EXPERIMENT = OUTPUT = CONTROL = COUNTS_FILE = REGION = MASKER_FILE = '' #files
EXPERIMENT_FORMAT=PK
OPEN_EXPERIMENT=False
DEBUG=False
DISCARD=0
LABEL = 'noname'
OUTPUT_FORMAT=PK
OPEN_OUTPUT=False
ROUNDING=False
CONTROL_FORMAT=None
REGION_FORMAT=BED
OPEN_REGION= False
FRAG_SIZE = 0
TAG_LENGTH = 0
SPAN=40
P_VALUE=0.01
HEIGHT_LIMIT=400
CORRECTION=1.
NO_SUBTRACT = False
DO_NORMALIZE = False
SPLIT_PROPORTION=0.1
SPLIT_ABSOLUTE=0
TRIM_PROPORTION=0.3
OPEN_CONTROL=False
NO_SORT=False
DUPLICATES=0
THRESHOLD=0
TRIM_ABSOLUTE=0
MAX_DELTA=250
MIN_DELTA=20
HEIGHT_FILTER=8
DELTA_STEP=1
VERBOSE=True
SPECIES='hg19'
CACHED=True
REPEATS=100
MAX_CORRELATIONS=200
KEEP_TEMP = False
POISSONTEST = 'height'
STRANDED_ANALYSIS = False
REMLABELS = ''
PROXIMITY=50
POSTSCRIPT=False
SHOWPLOTS=False
PLOT_PATH=None
LABEL1=""
LABEL2=""
BINSIZE=0.3
ZSCORE = 2
SDFOLD = 1
BLACKLIST=None
RECALCULATE=False
REGION_MINTAGS = 6
WINDOW_STEP = 0.1
POISSON_OPTIONS=("height", "numtags", "length")
TEMPDIR=[]


#Enrichment
PSEUDOCOUNT=False
LEN_NORM=False
TMM_NORM=False
N_NORM=False
SKIP_HEADER=False
TOTAL_READS_A=None
TOTAL_READS_B=None
TOTAL_READS_REPLICA=None
A_TRIM=0.05
M_TRIM=0.25
SKIP_PLOT=False
USE_REPLICA=False

#CONSTANTS
PLUS_STRAND = "+"
MINUS_STRAND = "-"
NO_STRAND = "."
EPSILON=1.0842021724855044e-19 #The smallest number above 0. Got from running 1./sys.maxint

NORMALIZE = 'normalize'
EXTEND = 'extend'
SUBTRACT = 'subtract'
SPLIT = 'split'
TRIM = 'trim'
FILTER = 'filter'
POISSON = 'poisson'
NOWRITE = 'nowrite'
DISCARD_ARTIFACTS = 'discard'
REMOVE_REGION = 'remove_regions'
REMOVE_DUPLICATES = 'remove_duplicates'
MODFDR = 'modfdr'
STRAND_CORRELATION = 'strand_correlation'

USE_MA = 'use_ma'
ENRICHMENT = 'enrichment'
CALCZSCORE = 'zscore'
PLOT='plot'






