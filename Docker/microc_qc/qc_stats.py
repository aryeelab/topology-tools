#!/usr/bin/env python3

'''
'''

from sys import argv
from matplotlib import pyplot as plt
import numpy
import argparse
from tabulate import tabulate



parser = argparse.ArgumentParser()
parser.add_argument('-p', help="Pairtools stat output")
parser.add_argument('-i', help="Pairtools Pairs file")
parser.add_argument('-d', help="Sample id")

args = parser.parse_args()


def qc_stats(stats_file):
	'''
	'''

	output_dict = {}
	with open(stats_file,'r') as f:
		for line in f:
			attrs = line.split()
			output_dict[attrs[0]] = attrs[1]

	table = []
	total_reads = int(output_dict["total"])
	total_reads_str = format(total_reads, ",d")
	table.append(["Total Read Pairs", total_reads_str, "100%"])

	unmapped_reads = int(output_dict["total_unmapped"])
	percent_unmapped = round(unmapped_reads*100.0/total_reads,2)
	unmapped_reads = format(unmapped_reads,",d")
	table.append(["Unmapped Read Pairs", unmapped_reads, f"{percent_unmapped}%"])

	mapped_reads = int(output_dict["total_mapped"])
	percent_mapped = round(mapped_reads*100.0/total_reads,2)
	mapped_reads_str = format(mapped_reads,",d")
	table.append(["Mapped Read Pairs", mapped_reads_str, f"{percent_mapped}%"])

	dup_reads = int(output_dict["total_dups"])
	percent_dups = round(dup_reads*100.0/total_reads,2)
	dup_reads = format(dup_reads,",d")
	table.append(["PCR Dup Read Pairs", dup_reads, f"{percent_dups}%"])

	nodup_reads = int(output_dict["total_nodups"])
	percent_nodups = round(nodup_reads*100.0/total_reads,2)
	nodup_reads_str = format(nodup_reads,",d")
	table.append(["No-Dup Read Pairs", nodup_reads_str, f"{percent_nodups}%"])

	cis_reads = int(output_dict["cis"])
	percent_cis = round(cis_reads*100.0/nodup_reads,2)
	cis_reads_str = format(cis_reads,",d")
	table.append(["No-Dup Cis Read Pairs", cis_reads_str, f"{percent_cis}%"])

	trans_reads = int(output_dict["trans"])
	percent_trans = round(trans_reads*100.0/nodup_reads,2)
	trans_reads = format(trans_reads,",d")
	table.append(["No-Dup Trans Read Pairs", trans_reads, f"{percent_trans}%"])

	cis_gt1kb = int(output_dict["cis_1kb+"])
	cis_lt1kb = cis_reads - cis_gt1kb
	percent_cis_lt1kb = round(cis_lt1kb*100.0/nodup_reads,2)
	percent_cis_gt1kb = round(cis_gt1kb*100.0/nodup_reads,2)
	cis_gt1kb = format(cis_gt1kb,",d")
	cis_lt1kb = format(cis_lt1kb,",d")
	valid_read_pairs = int(output_dict["trans"]) + int(output_dict["cis_1kb+"])
	percent_valid_read_pairs = round(valid_read_pairs*100.0/nodup_reads, 2)
	valid_read_pairs = format(valid_read_pairs,",d")
	table.append(["No-Dup Valid Read Pairs (cis >= 1kb + trans)", valid_read_pairs, f"{percent_valid_read_pairs}%"])
	table.append(["No-Dup Cis Read Pairs < 1kb", cis_lt1kb, f"{percent_cis_lt1kb}%"])
	table.append(["No-Dup Cis Read Pairs >= 1kb", cis_gt1kb, f"{percent_cis_gt1kb}%"])

	cis_gt10kb = int(output_dict["cis_10kb+"])
	percent_cis_gt10kb = round(cis_gt10kb*100.0/nodup_reads,2)
	cis_gt10kb = format(cis_gt10kb,",d")
	table.append(["No-Dup Cis Read Pairs >= 10kb", cis_gt10kb, f"{percent_cis_gt10kb}%"])

	cis_gt20kb = int(output_dict["cis_20kb+"])
	percent_cis_gt20kb = round(cis_gt20kb*100.0/nodup_reads,2)
	cis_gt20kb = format(cis_gt20kb,",d")
	table.append(["No-Dup Cis Read Pairs >= 20kb", cis_gt20kb, f"{percent_cis_gt20kb}%"])

	return table


def qc_histogram(pairfile):
	'''
	'''
	with open(pairfile, 'r') as f:
		total_pairs = 0; distances = []
		dist_20kb = 0; dist_10kb = 0
		trans_interactions = 0
		for line in f:
			if not line.startswith('#'):
				total_pairs+=1
				newline = line.rstrip('\n').split('\t')
				chr1 = newline[1];chr2 = newline[3]
				pos1 = newline[2];pos2 = newline[4]
				
				if chr1==chr2:
					distance = abs(int(pos2)-int(pos1))
					if distance == 0:
						distance = 1
					distances.append(distance)
					if distance >= 10000:
						dist_10kb += 1
					if distance >= 20000:
						dist_20kb += 1
				else:
					trans_interactions += 1

		# Plot the histogram
		n, bins, patches = plt.hist(numpy.log10(distances), 100, alpha=0.5)

		plt.xlabel('Log10 distances')
		plt.ylabel('Number of pairs')
		plt.title(r'Pair distances distribution')

# Save the histogram
		plt.savefig('hist.png')
		plt.show()

	return None


def create_html(sample_name, qc_table):
	'''
	'''
	page_title_text='MicroC QC Report: ' + sample_name
	title_text = 'MicroC QC Report: ' + sample_name
	text = 'This is a contact pairs microC report. Webpage for reference numbers: https://micro-c.readthedocs.io/en/latest/library_qc.html. \
	Need to change the values to represent percentages of the total number of reads (?)'
	hist_text = 'Histogram of valid contact pairs'

	# 2. Combine them together using a long f-string
	html = f'''
	    <html>
	        <head>
	            <title>{page_title_text}</title>
	        </head>
	        <body>
	            <h1>{title_text}</h1>
	            <p>{text}</p>
	            {tabulate(qc_table, tablefmt='html', stralign='center').replace('<table>',"<table border='1' cellpadding='4' cellspacing='7'>")}
	            <h2>{hist_text}</h2>
	            <img src='hist.png' width="700">
	        </body>
	    </html>
	    '''
	# 3. Write the html string as an HTML file
	with open('html_report.html', 'w') as f:
		f.write(html)
	return None




if __name__ == "__main__":

	qc_table = qc_stats(args.p)
	qc_histogram(args.i)
	create_html(args.d, qc_table)


