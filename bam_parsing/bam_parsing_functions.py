#!/usr/bin/env python
import pysam
import gzip
from io import StringIO
import sys

def extractUnmappedReads(input_bam,output_filepath):
    
    output_bam = pysam.AlignmentFile(output_filepath,"wb",template=input_bam)
    
    """
    Extract unmapped reads from a BAM file and write them into an output file.
    """
    print("Filtering unmapped reads:")
    for read in input_bam.fetch(until_eof=True):
        if read.is_unmapped:
            output_bam.write(read) 
            
    print("Closing connection")        
    output_bam.close()

    print("Indexing BAM file")
    pysam.index(output_filepath)
    
    
    input_bam.close()
    print("Done.")


def bamToFastQ(bam_file_path):
    """
    Convert BAM file to FastQ.gz reads file format
    """
    bam = pysam.AlignmentFile(bam_file_path,"rb")
    output_filename = bam_file_path.replace(".bam",".fastq.gz")
    
    # Redirect stdout to a temporary buffer
    original_stdout = sys.stdout
    sys.stdout = StringIO()
    
    # Run the function
    pysam.fastq(bam_file_path)
    
    # Get the content of the temporary buffer
    stdout_content = sys.stdout.getvalue()
    
    # Restore stdout
    sys.stdout.close()  # Close the StringIO object
    sys.stdout = original_stdout
    
    print("Converting BAM to FASTQ.gz")
    with gzip.open(output_filename,'wt') as output_fastq:
            #for read in bam.fetch(until_eof=True):
            output_fastq.write(stdout_content)
    
    output_fastq.close()            
print("Done.")
   
   

   
    
    