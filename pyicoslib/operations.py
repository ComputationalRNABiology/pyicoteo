Normalize = 'normalize'
Extend = 'extend'
Subtract = 'subtract'
Split = 'split'
Trim = 'trim'
Filter = 'filter'

class Poisson:
    incompatible_with = ['enrichment']
    
    def __str__(self):
        return 'poisson'

Convert = 'convert'
NoWrite = 'nowrite'
DiscardArtifacts = 'discard'
RemoveRegion = 'remove'
RemoveDuplicates = 'remove_duplicates'
ModFDR = 'modfdr'

class StrandCorrelation:
    def __str__(self):
        return 'strand_correlation'

class Enrichment:
    incompatible_with = ['modfdr', Filter, Poisson]
    def __str__(self):
        return 'enrichment'

class OperationFailed(Exception):
    pass
