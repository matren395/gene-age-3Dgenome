"""
Script to keep era genes, HARs, or UCEs that overlap 3D genomic features of interest.

"""

import argparse
import pybedtools as pbt

def get_args():
    # Catch command-line arguments
    parser = argparse.ArgumentParser(description="Script to keep era genes, HARs, UCEs that overlap 3D genomic features of interest")
    parser.add_argument("file", type=argparse.FileType('rU'), help='A file containing'
		'filenames for era genes, HARs, and UCEs you want to process,'
		'separated by newlines')
    parser.add_argument('-e', '--org', type=argparse.FileType('rU'), help='A file containing'
        'filenames for 3D genomic features you want to process,'
        'separated by newlines')
    return parser.parse_args()

    
    

def get_filenamesA(args):
	aFilenamesA = [line.strip() for line in args.file]
	return aFilenamesA
 
def get_filenamesB(args):
    aFilenamesB = [line.strip() for line in args.org]
    return aFilenamesB
		                         

def getFeatures(strFileName):
    btFeatures = pbt.BedTool(strFileName)
    return btFeatures

def sort(btObject):
    btSorted = btObject.sort()
    return btSorted

def saveBedTool(btObject, strFilename):
    btObject.saveas(strFilename)



def main():
    args = get_args()
    aFilenamesA = get_filenamesA(args)
    aFilenamesB = get_filenamesB(args)
    for filename in aFilenamesA:
        btfeatures = getFeatures(filename)
        for file in aFilenamesB:
            btOrg = getFeatures(file)
            btOverlap = btfeatures.intersect(btOrg, u=True)
            #sort
            btOverlap_Sorted = sort(btOverlap)
            #output saved
            saveBedTool(btOverlap_Sorted, '{0}_{1}_Overlap_Sorted.bed'.format(filename, file))
  

if __name__ == "__main__":
     main()
