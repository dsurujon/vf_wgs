
## DS 02/20/18
#  Filter population mutations by comparing against adaptation
#  experiments done in CDM (control for background mutations).
#  Retain only mutations that are present in one condition.
#
#  The experiment sheet provides a list of gd files to be analyzed.
#  It has two columns labeled "File" and "Group". The entries under
#  file should correspond to the .gd file names (.gd extension 
#  included), and Group could be either "E" (experimental) or "C"
#  (control).
#  The genbenk reference genome is used for adding locus tags to
#  mutations if they appear in a gene. 
#  
#  The filtering is done by keeping all mutations that satisfy the
#  following two conditions:
#    1. The control populations do not have >(lower cutoff) frequency
#    2. There is at least one experimental population that has
#       >(higher cutoff) frequency


import os
import glob
from Bio import SeqIO
import numpy as np
import pandas as pd
from optparse import OptionParser
import re
import traceback

# this is specific for this analysis on V. fisheri
# BreSeq was run on the Genbank reference genome, 
GB_to_RS = {'CP000020': 'NC_006840.2',
            'CP000021': 'NC_006841.2',
            'CP000022': 'NC_006842.1'}

# Use example:
# python filter_gd.py -i testdir -s Filter_Kan.csv -g NC_012469.gbk -o testout_10_50.csv -l 10 -u 50

options = OptionParser(usage='%prog -i inputdir -s exptsheet -g genome -o output -l lowercutoff -u uppercutoff [-C] [-M]',
                       description="Specify input directory, experiment sheet, reference genome file and output file, and cutoffs for frequency filtering. Optional: if you want to run this for clonal samples, add '-C' to the command. If you want to generate a table to use with the Muller_diagram package, add the -M option.")

options.add_option("-i","--inputdir",dest="inputdir",
                   help="input directory containing gd files")
options.add_option("-s","--infile",dest="inputfile",
                   help="input file (.csv) listing which gd files are experimental and which are control")
options.add_option("-g","--gbk",dest="genomefile",
                   help="reference genome file")
options.add_option("-l","--lcut",dest="lowercutoff",
                   help="lower cutoff for mutation frequency between [0,100]")
options.add_option("-u","--ucut",dest="uppercutoff",
                   help="higher cutoff for mutation frequency between [0,100]")
options.add_option("-o","--outfile",dest="outputfile",
                   help="output file (.csv)")
options.add_option("-C", action="store_true", dest="clonal")
options.add_option("-M", action="store_true", dest="makemullertable")



#read gd file, return a list of lines that have specific types of
#mutations (eg. "SNP")  
def reformat_gd(filename, reference, isclonal):
    filebasename = os.path.splitext(os.path.basename(filename))[0]
    f = open(filename)
    lines = [l.split() for l in f.readlines()]
    f.close()

    relevantlines = []
    relevantmutations = ["SNP","SUB","INS","DEL","MOB","AMP","CON","INV"]

    type=[]
    contig=[]
    pos=[]
    fromnt=[]
    tont=[]
    frequency=[]
	# for population (i.e. nonclonal samples) use the frequency listed on the gd file
    if isclonal==False:
        for line in lines:
            mut_type = line[0]
            if mut_type in relevantmutations:
                try:
                ## if mut_type is INS, need to get line[6]
                    if mut_type=="INS":
                        frequency.append(float(line[6].split('=')[1])*100)
                    else:
                        frequency.append(float(line[-1].split('=')[1])*100)
                except ValueError:
                    frequency.append(0)
                p=int(line[4])
                type.append(mut_type)
                pos.append(p)
                fromnt.append(reference[GB_to_RS[line[3]]].seq[p-1])
                tont.append(line[5])
                contig.append(line[3])
    # for clonal samples, add a frequency of 1 for each mutation
    else:
        for line in lines:
            mut_type = line[0]
            if mut_type in relevantmutations:
                frequency.append(100)
                p=int(line[4])
                type.append(mut_type)
                pos.append(p)
                fromnt.append(reference[GB_to_RS[line[3]]].seq[p-1])
                tont.append(line[5])
                contig.append(line[3])
    d = {'Contig': pd.Series(contig),
         'Position': pd.Series(pos),
         'From': pd.Series(fromnt),
         'To': pd.Series(tont),
         'Type': pd.Series(type),
         'Frequency': pd.Series(frequency)}
    df = pd.DataFrame(d)
    df=df[['Contig','From','Position','To','Type','Frequency']]
    df.rename(columns={'Frequency':'Frequency.'+filebasename}, inplace=True)
    return(df)

#merge dataframes from all population files
def make_merged_FT(filenames,reference, isclonal):
    print('reading', filenames[0])
    dfmerge=reformat_gd(filenames[0],reference,isclonal)
    numfiles=len(filenames)
    for i in range(1,numfiles):
        print('reading ', filenames[i])
        
        dfmerge=dfmerge.merge(reformat_gd(filenames[i],reference,isclonal),how="outer", on=["Contig","Position","From","To","Type"])
    return(dfmerge)



## check where each position lies. If there is a feature at that coordinate, get the locus tag
def get_locus_tags(mydf,feature_list):
    pos=[i for i in mydf["Position"]]
    lt=["" for i in pos]
    for i in range(0,len(pos)):
        for feature_dictionary in feature_list:
            p_i=pos[i]
            if feature_dictionary["start"]<=p_i and p_i<=feature_dictionary["end"]:
                lt[i]=feature_dictionary["locus_tag"][0]
                break
    mydf['LocusTag']=lt
    return(mydf)


#compare to control set (usually CDM)
#mydf: experiment dataframe
#cdmdf: control df
#fc_low/high: frequency cutoffs for control and experiment
#fl: featurelist
#filename: outputfilename
def make_ctrl_comparison(mydf,cdmdf,fc_low,fc_high,fl,filename):
    mergedf = cdmdf.merge(mydf, how="outer", on=['Contig','Type',"From","To","Position"],suffixes=[".CTRL",""])
    mergedf = mergedf.drop_duplicates()
    mergedf = mergedf.fillna(0)
    cdmcols = [col for col in mergedf.columns if 'CTRL' in col]
    exptcols = [col for col in mergedf.columns if 'EXPT' in col]
    df_sub=mergedf[(mergedf[cdmcols]<fc_low).all(axis=1) & (mergedf[exptcols]>fc_high).any(axis=1)]
    #df_sub=df_sub[(df_sub[exptcols]>fc_high).sum(axis=1)>1]
    df_multi=get_locus_tags(df_sub,fl)
    print(len(df_sub),len(set(df_multi['LocusTag'])))
    df_multi.to_csv(filename)
    return(df_multi)
    


def read_expt_sheet(filename):
    mylist = pd.read_csv(filename)
    exptgroup = list(mylist[mylist['Group']=="E"]['File'])
    ctrlgroup = list(mylist[mylist['Group']=="C"]['File'])
    return(exptgroup, ctrlgroup)

def main():
    # read input arguments
    opts, args = options.parse_args()
    genomefile = opts.genomefile
    inputdir = opts.inputdir
    inputfile = opts.inputfile
    outfile = opts.outputfile
    cut_low = int(opts.lowercutoff)
    cut_high = int(opts.uppercutoff)
    isclonal = opts.clonal
    makemullertable = opts.makemullertable
    if isclonal==None:
        isclonal=False
    if makemullertable==None:
        makemullertable=False
    
    # go to specified directory
    os.chdir(inputdir)
    # read reference genome
    # in order to support genbank with multiple records/contigs, read as dict
    strain_dict = SeqIO.to_dict(SeqIO.parse(genomefile, "genbank"))
    # read experiment sheet
    expts, ctrls = read_expt_sheet(inputfile)
    
    # Control data frame for mutation frequencies
    try:
        ctrlDF=make_merged_FT(ctrls,strain_dict, isclonal).drop_duplicates()
        ctrlDF.rename(columns=lambda x: x.replace('Frequency', 'CTRL'), inplace=True)
        # Experiment data frame for mutation frequencies
        exptDF=make_merged_FT(expts,strain_dict,isclonal).drop_duplicates()
        exptDF.rename(columns=lambda x: x.replace('Frequency', 'EXPT'), inplace=True)
        print('Done reading gd files')
    except IndexError:
        if isclonal==False:
            print('Running in population mode. Check if your samples are clonal. If so, please re-run the script with the additional tag -C')
            print(traceback.format_exc())
            print('Quitting...')
            quit()

    ## generate a feature dictionary for the reference genome
    feature_list=[]
    for record_id, record_seq in strain_dict.items():
        for feature in record_seq.features:
            if feature.type=="CDS":
                fd={"locus_tag":feature.qualifiers["locus_tag"],
                    "start": int(feature.location.start),
                    "end": int(feature.location.end)}
                feature_list.append(fd)
    print('Done making feature lists')

    print("High cutoff: %s, Low cutoff:%s\n" %(cut_low, cut_high))
    filtered_df = make_ctrl_comparison(exptDF,ctrlDF,cut_low,cut_high,feature_list,outfile)
    
    if makemullertable:
        filtered_df['Population']='X'
        filtered_df['Chromosome']=1
        filtered_df['Trajectory']=[i for i in range(filtered_df.shape[0])]
        filtered_df.rename(columns={'Type':'Class'}, inplace=True)
        filtered_df['Mutation'] = filtered_df['From']+'>'+filtered_df['To']
        exptcolnames = ['EXPT.'+os.path.splitext(os.path.basename(filename))[0] for filename in expts]
        relevantcols = ['Population', 'Trajectory', 'Chromosome', 'Position', 'Class', 'Mutation']+exptcolnames
        muller_df = filtered_df[relevantcols]
        exptrename = {exptcol:re.search(r'\d+', exptcol).group() for exptcol in exptcolnames}
        muller_df.rename(columns = exptrename, inplace=True)
        muller_filename = os.path.splitext(outfile)[0]+'_Muller.csv'
        muller_df.to_csv(muller_filename, index=False)

if __name__ == '__main__':
    main()
