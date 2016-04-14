from io import *
from glob import *
import numpy

def collapseWeightsToCSV(folder:str, outfile:str):
	output = open(outfile, 'w');

	files = glob(folder + "*.txt");
	for file in files:
		rawfile = open(file, 'r');
		for line in rawfile:
			entries = line.split();
			if entries[0] == 'Weights:':
				w = numpy.zeros([len(entries) - 1]);
				for i in range(1, len(entries)):
					w[i - 1] = float(entries[i]);
				w = w / numpy.sum(w);
				for i in range(w.shape[0]):
					output.write(str(w[i]));
					output.write(',');
				output.write('\n');
				break;
		rawfile.close();
	output.close();

def collapseParamsToCSV(folder:str, outfile:str):
	output = open(outfile, 'w');
	output.write("Ensemble Size,Backrub,Boltzmann,Steepness\n");

	files = glob(folder + "*.txt");
	for file in files:
		rawfile = open(file, 'r');
		for i in range(4):
			line = rawfile.readline();
			entries = line.split();
			output.write(entries[-1].replace('\n', '') + ',');
		output.write('\n');
		rawfile.close();
	output.close();