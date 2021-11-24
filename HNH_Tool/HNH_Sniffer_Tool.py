'''
Tool for finding HNH motifs in genomes, particularly in bacteriophage.
Caleb A. Collar
BIO_400Q
Notes:
- Requres BLAST+ local installation.
- No databases needed.
- Works with Python 3x
'''
#Imports
from Bio.Blast.Applications import NcbiblastpCommandline
import io
from Bio.Blast import NCBIXML
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
#Basic UI to get user's file for comparison.
userIn = input("Please input name of file within the FASTA_DIR folder to check 'ex. MyPhage.faa':")
# Run BLAST and output to xml
output = NcbiblastpCommandline(query="HNH_Motif/HNH2.faa", subject="FASTA_DIR/" + userIn, outfmt=5)()[0]
print("Working") #Works until here to parse results in xml
print('\n')
result_handle = NCBIXML.read(io.StringIO(output))
# Print information on the result
E_VALUE_THRESH = 0.04 #e.x. 'e = 0.05' is 5 in 100 chance of being random.
num_hits = 0
hnh_hits = 0
#terminase_hits = 0
for alignment in result_handle.alignments:
    for hsp in alignment.hsps:
        num_hits = num_hits + 1
        if hsp.expect < E_VALUE_THRESH:
            hnh_hits = hnh_hits + 1
            print('------------Possible Hit With HNH Motif------------')
            print(alignment.title)
            print('e-Value:', hsp.expect)
            print('Length:', alignment.length)
            print('Position in current gene:', hsp.sbjct_start, '-', hsp.query_end)
            #print('Identities:', hsp.identities)
            #print('Positives:', hsp.positives)
            #print('Gaps:', hsp.gaps)
            #print('Bits:', hsp.bits)
            print('\n')
            print('Visualization of alignment:')
            print(hsp.query[0:75] + '...')
            print(hsp.match[0:75] + '...')
            print(hsp.sbjct[0:75] + '...')
            print('---------------------------------------------------')
            print('\n')
        else:
            if num_hits <= 1:
                print('No hits over e threshold of', E_VALUE_THRESH, 'from file:', 'FASTA_DIR/' + userIn)
#End Info
if hnh_hits > 0:
    print('Since the tool found a hit with an e-Value high enough, it is likely that this sequence contains an HNH motif.')
print('Done.')
print('\n')