## DS 08/22/19
# updated 12/17/24
# Copy raw BreSeq output gd files into a new directory
# and rename them from "output.gd" to "[sampleame].gd" 
# where [samplename] is the name of the parent directory 
# of the output.gd file under ExpOut

import os
import shutil
from optparse import OptionParser

## Use example 
# python copy_gd.py -e expID -o outputdir

options = OptionParser(usage='%prog -e expID -o outputdir -t filetype',
		description = "Specify the local raw output directory (where the gd/vcf files are currently) and the directory you want to copy them to. The output.gd or .vcf files will automatically be renamed in the new directory as their sample name.")
		
options.add_option("-e", "--expID", dest="expID", help="local directory with raw breseq output")
options.add_option("-o", "--outputdir", dest="outputdir", help="output directory")
options.add_option("-t", "--filetype", dest="filetype", help="output filetype (vcf or gd)")


def main():
	#read input args
	opts, args = options.parse_args()
	expID_dir = opts.expID
	outdir = opts.outputdir
	file_ext = opts.filetype
	
	# if output directory doesn't exist, make one
	if os.path.exists(outdir)==False:
		os.makedirs(outdir)
	
	# list all subdirectories under expID/Out
	raw_output_dir = expID_dir
	samples = os.listdir(raw_output_dir)
	
	# for each sample, copy the gd file to outdir, and rename it
	for samplename in samples:
		if samplename!="Bams":
			print("Transferring" , samplename)
			try:
				source_gd = os.path.join(raw_output_dir, samplename, "output", f'output.{file_ext}')
				shutil.copy(source_gd, outdir)
				oldgdname = os.path.join(outdir, f'output.{file_ext}')
				newgdname = os.path.join(outdir, samplename + "." + file_ext)
				os.rename(oldgdname,newgdname)
			except IOError:
				print("Error transferring file", samplename)

if __name__ == '__main__':
	main()
