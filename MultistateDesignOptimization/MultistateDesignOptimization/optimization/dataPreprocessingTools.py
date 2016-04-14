from io import *
import numpy

def targetFASTAstrip(source, indices, out):
	"""
	Reads SOURCE, writes only indices defined by INDICES to OUT.
	Each index is the first entry in its own line.
	"""
	infile = open(source, 'r');
	index = open(indices, 'r');
	outfile = open(out, 'w');

	positions = [];
	for line in index:
		entries = line.split(' ');
		positions.append(int(entries[0]));

	for line in infile:
		if line[0] == '>':
			outfile.write(">NULL\n");
		else:
			for i in range(len(positions)):
				outfile.write(line[positions[i]]);
			outfile.write('\n');

	infile.close();
	index.close();
	outfile.close();

def extractModelPositions(source, positions, out):
	"""
	Extracts the set of positions defined in POSITIONS from a .dat file
	"""
	
	infile = open(source, 'r');
	outfile = open(out, 'w');

	firstLine = True;
	lowPos = 0;
	newPos = numpy.zeros([len(positions)]);
	for i in range(len(positions)):
		newPos[i] = i;

	for line in infile:
		if firstLine:
			firstLine = False;
			entries = line.split(' ');
			lowPos = int(entries[0]);
			# new lowest position is 0
			outfile.write("0 ");
			outfile.write("{:d} ".format(len(positions)));
			outfile.writelines(entries[2]);
			continue;

		entries = line.split(' ');
		if int(entries[2]) - lowPos in positions:
			i = 0;
			for entry in entries:
				if i != 2:
					outfile.write(entry);
					if entry[-1] != '\n':
						outfile.write(' ');
				else:
					outfile.write("{:d} ".format(positions.index(int(entries[2]) - lowPos)));
				i += 1;

	infile.close();
	outfile.close();

	return None

def extractFASTAPositions(source, positions, out):
	"""
	Extracts a set of positions, defined by the list POSITIONS, from a FASTA
	file.
	"""

	infile = open(source, 'r');
	outfile = open(out, 'w');

	for line in infile:
		if line[0] == '>':
			outfile.write(">Null\n");
		else:
			for i in range(len(positions)):
				outfile.write(line[positions[i]]);
			outfile.write("\n");

	infile.close();
	outfile.close();

	return None

def glycineDeweighter(source:str, amount = 0.75):

	name = source.split('.');
	name = name[0] + " G deweighted.dat";

	infile = open(source, 'r');
	outfile = open(name, 'w');

	firstLine = True;

	for line in infile:

		if firstLine:
			firstLine = False;
			outfile.write(line);
			continue;

		entries = line.split(' ');
		for i in range(22):
			if i != 8:
				outfile.write(entries[i]);
			else:
				outfile.write("{:.6f}".format(amount * float(entries[i])));
			outfile.write(' ');
		outfile.write(entries[22]);

	infile.close();
	outfile.close();
