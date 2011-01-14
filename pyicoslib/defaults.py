#Formats
ELAND = 'eland'
ELAND_EXPORT = 'eland_export'
BED = 'bed'
WIG = 'bed_wig'
VARIABLE_WIG = 'variable_wig'
FIXED_WIG = 'fixed_wig'
PK = 'bedpk'
SPK = 'bedspk'
SAM = 'sam'
CLUSTER_FORMATS = (WIG, VARIABLE_WIG, FIXED_WIG, PK, SPK)
WIG_FORMATS = (WIG, VARIABLE_WIG, FIXED_WIG)
READ_FORMATS = (ELAND, BED, WIG, PK, SPK, SAM) #formats that we actually can read as

PLUS_STRAND = "+"
MINUS_STRAND = "-"
NO_STRAND = "."

#Default values for parser
EXPERIMENT=''
EXPERIMENT_FORMAT=PK
OPEN_EXPERIMENT=False
DEBUG=False
DISCARD=0
OUTPUT=''
CONTROL=''
LABEL = 'noname'
OUTPUT_FORMAT=PK
OPEN_OUTPUT=False
ROUNDING=False
CONTROL_FORMAT=None
REGION=''
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
SPLIT_PROPORTION=0.9
SPLIT_ABSOLUTE=0
TRIM_PROPORTION=0.3
OPEN_CONTROL=False
NO_SORT=False
DUPLICATES=3
THRESHOLD=0
TRIM_ABSOLUTE=0
MAX_DELTA=500
MIN_DELTA=20
HEIGHT_FILTER=8
DELTA_STEP=1
VERBOSE=True
SPECIES='hg19'
CACHED=True
REPEATS=100
MASKER_FILE=''
MAX_CORRELATIONS=200
KEEP_TEMP = False
POISSONTEST = 'height'
STRANDED_ANALYSIS = False
WRITE_FORMATS = (ELAND, BED, WIG, VARIABLE_WIG, PK, SPK) #formats we can actually write as
REMLABELS = ''
PROXIMITY=50
POSTSCRIPT=False
SHOWPLOTS=False

"""The minimum number of overlapping positive and negative strand reads to include them in the correlation calculation"""

NORMALIZE = 'NORMALIZE'
EXTEND = 'extend'
SUBTRACT = 'subtract'
SPLIT = 'split'
TRIM = 'trim'
FILTER = 'filter'
POISSON = 'poisson'
NOWRITE = 'nowrite'
DISCARD_ARTIFACTS = 'discard'
REMOVE_REGION = 'remove'
REMOVE_DUPLICATES = 'remove_duplicates'
ModFDR = 'modfdr'
STRAND_CORRELATION = 'strand_correlation'
ENRICHMENT = 'enrichment'







