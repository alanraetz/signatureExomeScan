
from filesplit.split import Split

split = Split("exomeScanPeptides.fasta", "split")

split.bylinecount(20000) 