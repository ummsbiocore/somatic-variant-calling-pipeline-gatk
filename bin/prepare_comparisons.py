#!/usr/bin/env python3

import sys, os, argparse, csv, shutil
from collections import defaultdict

def main(args):

	if os.path.basename(args.groups)[:7] != 'NO_FILE' and os.path.basename(args.comparisons)[:7] != 'NO_FILE': 

		validate_files(args)
	
		comparisons = read_comparisons(args.comparisons)
		
		for comparison in comparisons:
			process_comparison(comparison)

def validate_files(args):

	warnings = []
	errors = []

	warnings.extend(check_filenames([args.groups, args.comparisons]))
	comparisons = check_comparison_file(args.comparisons, errors)
	samples, groups = check_group_file(args.groups, errors)

	for comparison in comparisons:
		treat, control, name, column = comparison
		
		if column not in groups.keys():
			errors.append('\tError: %s from comparison file not a column in groups file.' % (column))
		if treat not in groups[column]:
			errors.append('\tError: %s from comparison file not in %s from groups file.' % (treat, column))
		if control not in groups[column]:
			errors.append('\tError: %s from comparison file not in %s from groups file.' % (control, column))

	if len(warnings) > 0:
		print('%d Warnings Detected:\n' % len(warnings))
		print('\n\n'.join(warnings))

	if len(errors) > 0:
		error_message = '%d Errors Detected in Input Files:\n\n%s' % (len(errors), '\n\n'.join(errors))
		sys.exit(error_message)

def check_filenames(input_files):
	warnings = []

	problem_files = set()
	for input_file in input_files:
		if " " in input_file:
			problem_files.add(input_file)		
	
	if len(problem_files) > 0:
		warnings.append("\tWarning: Although allowed in the DE module, spaces are generally not recommended in filenames. Affected files: %s" % ('\n'.join((['\t\t%s' % input_file for input_file in problem_files]))))
		
	return(warnings)

def check_group_file(group_infile, errors):

	with open(group_infile) as infile:
		dialect = csv.Sniffer().sniff(infile.readline(), delimiters=',\t')
		infile.seek(0)
		reader = csv.reader(infile, dialect)
		header = next(reader)
		if len(header) <= 1:
			error_message = "Error: Only one tab-delimited or comma-delimited column detected in %s. Is this file either tab-delimited or comma-delimited?" % (group_infile)
			sys.exit(error_message)

		header = [i.replace('\ufeff', '').lower() for i in header]

		if 'sample_name' not in header:
			sys.exit('Error: "sample_name" must be a column in %s\n%s' % (group_infile, '\n'.join(errors)))

		samples = set()
		groups = defaultdict(set)

		column_key = {}
		for i, column in enumerate(header):
			column_key[i] = column
		num_columns = len(header)

		for line in reader:
			samples.add(line[0])
			for i, column in enumerate(line):
				groups[column_key[i]].add(column)

	return samples, groups

def check_comparison_file(comparison_infile, errors):

	comparisons = []
	
	with open(comparison_infile) as infile:
		dialect = csv.Sniffer().sniff(infile.readline(), delimiters=',\t')
		infile.seek(0)
		reader = csv.reader(infile, dialect)
		header = next(reader)
		
		if len(header) <= 1:
			error_message = "Error: Only one tab-delimited or comma-delimited column detected in %s. Is this file either tab-delimited or comma-delimited?" % (comparison_infile)
			sys.exit(error_message)

		header = [i.replace('\ufeff', '').lower() for i in header]

		problem = False

		if 'controls' in header:
			controls_column = header.index('controls')
		else:
			errors.append('\tError: Cannot find "controls" column in %s' % comparison_infile)
			problem = True

		if 'treats' in header:
			treats_column = header.index('treats')
		else:
			errors.append('\tError: Cannot find "treats" column in %s' % comparison_infile)
			problem = True

		if 'names' in header:
			names_column = header.index('names')
		else:
			errors.append('\tError: Cannot find "names" column in %s' % comparison_infile)
			problem = True

		if problem:
			sys.exit('\n'.join(errors))

		if 'grouping_column' in header:
			grouping_column = header.index('grouping_column')

		for line in reader:

			if len(line) == 3:
				comparisons.append((line[treats_column], line[controls_column], line[names_column], 'group'))
			elif len(line) == 4:

				if 'grouping_column' not in header:
					sys.exit('\tError: When using four column approach, "grouping_column" must exist in %s' % (comparison_infile))

				if line[grouping_column] == '':
					comparisons.append((line[treats_column], line[controls_column], line[names_column], 'group'))
				else:
					comparisons.append((line[treats_column], line[controls_column], line[names_column], line[grouping_column]))
			else:
				error_message = '\tError: Incorrect number of columns in comparison file. 3 or 4 required. %d columns in this line: \n\t\t%s' % (len(line), '\n'.join(line))
				sys.exit(error_message)
	
	return(comparisons)

def read_comparisons(comparison_infile):

	comparisons = []

	with open(comparison_infile) as infile:
		dialect = csv.Sniffer().sniff(infile.readline(), delimiters=',\t')
		infile.seek(0)
		reader = csv.reader(infile, dialect)
		header = next(reader)
		

		header = [i.replace('\ufeff', '').lower() for i in header]


		control_column = header.index('controls')
		treats_column = header.index('treats')
		names_column = header.index('names')

		for line in reader:

			if len(line) == 3:
				comparisons.append((line[treats_column], line[control_column], line[names_column], 'group'))
			else:
				grouping_column = header.index('grouping_column')
				if line[grouping_column] == '':
					comparisons.append((line[treats_column], line[control_column], line[names_column], 'group'))
				else:
					comparisons.append((line[treats_column], line[control_column], line[names_column], line[grouping_column]))
		
		return comparisons

def process_comparison(comparison):
	
	treatment_group, control_group, name, column = comparison

	groups = read_groups(args.groups, column)

	os.makedirs('%s/controls/' % name, exist_ok=True)
	os.makedirs('%s/treatments/' % name, exist_ok=True)

	for f in groups[control_group]:
		shutil.copy('%s_recal.bam' % (f), '%s/controls/' % (name))
		shutil.copy('%s_recal.bam.bai' % (f), '%s/controls/' % (name))

	for f in groups[treatment_group]:
		shutil.copy('%s_recal.bam' % (f), '%s/treatments/' % (name))
		shutil.copy('%s_recal.bam.bai' % (f), '%s/treatments/' % (name))
		

def read_groups(group_infile, column):

	groups = defaultdict(list)

	with open(group_infile) as infile:
		dialect = csv.Sniffer().sniff(infile.readline(), delimiters=',\t')
		infile.seek(0)
		reader = csv.reader(infile, dialect)
		header = next(reader)

		header = [i.replace('\ufeff', '').lower() for i in header]

		sample_name_column = header.index('sample_name')

		if column.lower() in header:
			index = header.index(column)
		else:
			index = 1

		for line in reader:
			groups[line[index]].append(line[sample_name_column])

	return groups

def parseArguments():
	parser = argparse.ArgumentParser(prog="prepare_comparisons", description='', usage='%(prog)s [options]')
	
	input_args = parser.add_argument_group('Input')
	input_args.add_argument('-g', '--groups', required=True, help='Name of groups file.', metavar='', dest='groups')
	input_args.add_argument('-x', '--comparisons', required=True, help='Name of comparisons file.', metavar='', dest='comparisons')
	return parser.parse_args()

if __name__ == "__main__":
	args = parseArguments()
	main(args)