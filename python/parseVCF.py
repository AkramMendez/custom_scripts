#!/usr/bin/env python
import sys
import os
import argparse
import csv
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages



parser = argparse.ArgumentParser(description='Script for parsing a VCF file to generate a summary report and a plot the distribution of variants.')

# Set arguments for input and output files
parser.add_argument('input_file', type=str, help='Path to the input VCF file.')
parser.add_argument('output_path', type=str, default= os.getcwd(), help='Path for saving the results (Default: Current directory). A report (.TXT) and a plot (.PDF) are generated and saved into the defined output path.')

# Parse the command-line arguments
args = parser.parse_args()

vcf_path = os.path.abspath(args.input_file)
input_filename, _ = os.path.splitext(os.path.basename(vcf_path))
output_path = args.output_path
outpath_txt = ''.join([os.path.abspath(output_path),'/',input_filename,"_summary_report.txt"])
outpath_pdf = ''.join([os.path.abspath(output_path),'/',input_filename,"_summary_report.pdf"])

#------------------------------------


#------------------------------------

def parseLine(line,sep):
    if sep == '=':    
        metinfo = line[2:].split(sep)[1]
    else:
        metinfo = line[2:].split(sep)[1]

    return metinfo

def getMetadata(vcf_path):
    
    metadata = []
    header = []
    with open(vcf_path,'r') as vcf:
        n_lines = 0
        for line in vcf:
            if line.startswith('##'):
                metadata.append(line.strip())
            elif line.startswith('#'):
                header.append(line.strip('\t'))
            n_lines += 1
            
    number_records = n_lines - len(metadata) - len(header)

    # Extract metadata information for the short tags:
    tags = ['fileformat','fileDate','source','reference']
    metadata_info = []
    for tag in tags:
        for line in metadata:
            if tag in line:
                parsed_tag = parseLine(line,'=')      
        
        metadata_info.append(parsed_tag)
    
    # Extract metadata information for the longer tags:
    filter_line = [line for line in metadata if line.startswith("##FILTER")]
    info_lines = [line for line in metadata if line.startswith("##INFO")]
    format_line = [line for line in metadata if line.startswith("##FORMAT")]
    contig_line = [line for line in metadata if line.startswith("##contig")]
    
    # Get sample IDs
    samples = header[0].strip().split('\t')[9:]
    number_samples = len(samples)
    format = parseLine(format_line[0],'=<ID=').split(',')[0]
    filter = parseLine(filter_line[0],'=<ID=').split(',')[0]
    info = [parseLine(x,'=<ID=').split(',')[0] for x in info_lines]
    contig = parseLine(contig_line[0],'=<ID=').split('=')[0].replace(">","")

    for x in [number_samples,number_records,filter,contig,format,info]:
        metadata_info.append(x)
        
    return(metadata_info)

def getData(vcf_path):
    with open(vcf_path,'r') as vcf:
        data = []
        split_data = []
        for line in vcf:
            if not line.startswith('##'):
                data.append(line.strip())
                
    split_data = [line.split('\t') for line in data]
    split_data_df = pd.DataFrame(split_data[1:], columns=split_data[0])
    return(split_data_df)

vcf_data = getData(vcf_path)                
#vcf_data

sample_columns = vcf_data.columns[9:]
#sample_columns

def cleanRow(row):
    return row.replace(r'[a-zA-Z]','')

vcf_data
tmp  = vcf_data['INFO'].str.split(';', expand=True)
tmp.columns=['AC','AN','VT']
tmp['VT'] = tmp['VT'].str.replace('VT=','')

variant_types_count_df = tmp.groupby('VT').size().reset_index()
variant_types_count_df.columns = ['Variant Type','Count']
variant_types_count_df


transitions = {
        'A': 'G',
        'G': 'A',
        'C': 'T',
        'T': 'C'
}

transversions = {
        'A': ['C','T'],
        'G': ['C','T'],
        'C': ['A','G'],
        'T': ['A','G'],
}

def is_homozygous(samples):
        n_homozygous = samples.isin(['0|0','1|1']).sum()
        return n_homozygous

def is_heterozygous(samples):
        n_heterozygous = samples.isin(['0|1','1|0']).sum()
        return n_heterozygous

def is_transition(ref,alt,info):
        if len(ref)==1 and 'SNP' in info:
                return alt in transitions.get(ref)
        else:
                return False
        
def is_transversion(ref,alt,info):
        if len(ref)==1 and 'SNP' in info:
                return alt in transversions.get(ref)
        else:
                return False

def substitution_type(ref,alt,info):
        if 'SNP' in info and is_transition(ref,alt,info):
                return ''.join([ref,'>',alt])
        elif 'SNP' in info and is_transversion(ref,alt,info):
                return ''.join([ref,'>',alt])
        else:
                return False
        
def deletion_type(ref,alt,info):
        if 'INDEL' in info and len(ref) < len(alt):
                return 'Insertion'
        elif 'INDEL' in info and len(ref) > len(alt):
                return 'Deletion'
        else:
                return False

tmp1 = vcf_data
tmp1['homozygous'] = tmp1.apply(lambda x: is_homozygous(x[sample_columns]), axis=1)
tmp1['heterozygous'] = tmp1.apply(lambda x: is_heterozygous(x[sample_columns]), axis=1)
tmp1['transition'] = tmp1.apply( lambda x: is_transition(x['REF'],x['ALT'],x['INFO']), axis =1)
tmp1['transversion'] = tmp1.apply( lambda x: is_transversion(x['REF'],x['ALT'],x['INFO']), axis =1)
tmp1['substitution_type'] = tmp1.apply( lambda x: substitution_type(x['REF'],x['ALT'],x['INFO']), axis =1)
tmp1['indel_type'] = tmp1.apply( lambda x: deletion_type(x['REF'],x['ALT'],x['INFO']), axis =1)
#tmp1


metrics = pd.DataFrame(index=[0])
metrics['SNPs count:'] = variant_types_count_df['Count'][1]
metrics['INDELs count:'] = variant_types_count_df['Count'][0]
metrics['Homozygous alleles count:'] = tmp1['homozygous'].sum()
metrics['Heterozygous alleles count:'] = tmp1['heterozygous'].sum()
metrics['Transitions (Ti) count:'] = tmp1[tmp1['INFO'].str.contains('SNP') & tmp1['transition']].shape[0]
metrics['Transversions (Tv) count:'] = tmp1[tmp1['INFO'].str.contains('SNP') & tmp1['transversion']].shape[0]
metrics['Ti/Tv ratio:'] = round(metrics['Transitions (Ti) count:']/metrics['Transversions (Tv) count:'],2)
metrics['Insertions count:'] = tmp1[tmp1['INFO'].str.contains('INDEL') & tmp1['indel_type'].str.match('Insertion')].shape[0]
metrics['Deletions count:'] = tmp1[tmp1['INFO'].str.contains('INDEL') & tmp1['indel_type'].str.match('Deletion')].shape[0]
metrics['Indel ratio:'] = round(metrics['Insertions count:']/metrics['Deletions count:'],2)
metrics.reset_index(drop=True,inplace=True)
metrics.fillna(0)

# Count substitution types
substitution_type_counts = tmp1[tmp1['INFO'].str.contains('SNP')].groupby('substitution_type').size().sort_values(ascending=False).to_frame()
substitution_type_counts.reset_index(drop=False,inplace=True)
substitution_type_counts.columns = ['Substitution type','Count']
substitution_type_counts

# Analyse possible mutational processes that might be operative:

def mutProcess(metrics,substitution_type_counts):
    output = []
    output.append("# Possible Mutational Processes:")
    output.append("========================================")
    ts_tv_ratio = metrics['Ti/Tv ratio:'][0]
    
    if ts_tv_ratio is not None:
        if ts_tv_ratio > 0.5:
            output.append("Ti/Tv > 0.5: Possible bias towards transitions.")
        elif ts_tv_ratio < 0.5:
            output.append("Ti/Tv < 0.5: Possible bias towards transversions.")
        elif ts_tv_ratio == 0.5:
            output.append("Ti/Tv == 0.5: Expected Ti and Tv's distribution by chance.")
    else:
        output.append("Ti/Tv ratio not available in metrics.")
    
    
    a = substitution_type_counts.iloc[0][0]
    if transitions.get(a[0])==a[2]:
        output.append(f"The top most frequent variant: {a} is a Transition")
        if(a[0]=='C'):
            output.append(f"    Change '{a}' possibly due to deamination")
    elif a[2] in transversions.get(a[0]):
        output.append(f"The top most frequent variant: {a} is a Transversion")
        
    
    indel_ratio = metrics['Indel ratio:'][0]
    if indel_ratio > 1:
        output.append("Indel ratio > 1: More Insertions detected.")
    elif indel_ratio < 1:
        output.append("Indel ratio < 1: More Deletions detected.")
    elif indel_ratio ==1:
        output.append("Indel ratio == 1: Equal proportion of Ins/Del detected")
    
    return output

# Write report:
def decor_line():
    return report.write("========================================\n")

metadata_df  = getMetadata(vcf_path)

colnames = ['File Format:','File Date:','Source:','Reference:','Number of Samples:','Number of Records:','Filter:','Contig:','Format:','Info:']

#------------------------------------------------
# Saving metrics TXT report:
#------------------------------------------------
with open(outpath_txt,'w',newline='') as report:
    print("========================================")
    print(f"Saving TXT report to {outpath_txt}")
    
    report_writer = csv.writer(report, delimiter='\t')
    decor_line()
    report.write("#       Metadata Information              \n")
    decor_line()
    for colname,data in zip(colnames,metadata_df):
        report_writer.writerow([colname,data])
    decor_line()
    report.write("#  Summary Metrics:         \n")
    report.write("#  SNPs, INDELs, Homozygous, Heterozygous \n")
    decor_line()
        
    for _,data in metrics.transpose().iterrows() :
        report_writer.writerow([_,data[0]])
        #report_writer.writerow([_,data.iloc[0]])
    decor_line()
    report.write("#       Substitution types            \n")
    decor_line()
    for _,data in substitution_type_counts.iterrows():
        report_writer.writerow(list(data))
    
    decor_line()
    
    mutational_processes_report = pd.DataFrame(mutProcess(metrics,substitution_type_counts))
    for _,data in mutational_processes_report.iterrows():
        report_writer.writerow(data)
        
    decor_line()
    print(f"Report saved.")    

# Distribution of variants across the chromosome and possible mutational processes that might be operative.

# Plot SNPs and INDELs density:

def labelVariantType(df):
    if df['transition'] and 'SNP' in df['INFO']:
        return 'SNP (Ti)'
    elif df['transversion'] and 'SNP' in df['INFO']:
        return 'SNP (Tv)'
    elif 'Insertion' in df['indel_type'] and 'INDEL' in df['INFO']:
        return 'INDEL (Ins)'
    elif 'Deletion' in df['indel_type'] and 'INDEL' in df['INFO']:
        return 'INDEL (Del)'
    else:
        return False
        
    
df =  vcf_data[['POS','INFO','transition','transversion','indel_type']].copy()
df['VariantType'] = df.apply( lambda x: labelVariantType(x), axis =1)

bins = np.arange(0,int(df['POS'].max()),100000)

pos = df['POS'].values.astype(int)

def plotVariantDist(pos,variant_type,color, title, ax):
    
    window_size=1000000

    # split chromosome positions in 1 Mb bins:
    bins = np.arange(0, pos.max(), window_size)

    # use window midpoints as x coordinate
    x = bins[:-1]
    
    # compute variant density in each window
    count, _ = np.histogram(pos, bins=bins)

    y = count

    # plot
    sns.despine(ax=ax, offset=10)
    ax.plot(x, y, label=variant_type,color=color)
    ax.legend(loc='upper right')
    ax.set_xlabel('Chromosome position (Mb)')
    ax.set_ylabel('Number of Variants')
    if title:
        ax.set_title(title)
    else:
        ax.set_title('Variant Distribution')



variant_types = df['VariantType'].unique()
colors = sns.color_palette()

#------------------------------------------------
# Saving Variant Distribution Plot:
#------------------------------------------------
with PdfPages(outpath_pdf) as pdf:
    print("========================================")
    fig, ax = plt.subplots(figsize=(12, 3))
    fig.subplots_adjust(left=0.15, right=0.85, top=0.85, bottom=0.2)
    for i,variant in enumerate(variant_types):
        
        pos_i = df.loc[df['VariantType']==variant,'POS'].values.astype(int)
        color = colors[i]
        #plt.figure()
        plotVariantDist(pos=pos_i,variant_type=variant,color=color,title='',ax=ax)
    
    print(f"Saving Variant Distribution Plot to {outpath_pdf}")
    pdf.savefig()
    
print("========================================")
print("Done.")


    
