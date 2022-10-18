#!/usr/bin/env python3

'''
'''

from __future__ import division
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i', help="Pairtools Pairs file")
parser.add_argument('-r', help="Region")
parser.add_argument('-o', help="Output")

args = parser.parse_args()

def region_coverage(pairs_file, region, out_file):
	'''
	Parse the pairs file, count the number of pairs in the region
	of interest and output them in a new file.
	'''
	output = open(out_file, 'w')
	output.write('#Region selected: ' + region + '\n')
	output.write('#columns: readID chrom1 pos1 chrom2 pos2 strand1 strand2 pair_type\n')
	region=region.split()
	total_lines = 0; region_lines = 0
	with open(pairs_file, 'r') as file:
		for line in file:
			if not line.startswith('#'):
				total_lines += 1
				pair_data = line.split()
				if (pair_data[1] == region[0]) and (pair_data[3] == region[0]):
					if (int(pair_data[2]) >= int(region[1])) and (int(pair_data[2]) <= int(region[2])):
						if (int(pair_data[4]) >= int(region[1])) and (int(pair_data[4]) <= int(region[2])):
							region_lines += 1
							output.write(line)
	output.close()
	return total_lines, region_lines





if __name__ == "__main__":
	
	total_lines, region_lines = region_coverage(args.i, args.r, args.o)

	print(str(region_lines/total_lines))