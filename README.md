# DeepWHISTLE
prerequisite packages: tensorflow version 1.14.0, numpy, pandas, argparse 

inputs for the model: sequence no less than 1001bp long in fasta format, and genetic information from WHISTLE (http://180.208.58.19/whistle/index.html). all bases except the first and last 500bp will be tested to see if it is a prediction site

run python deepwhistle_restore.py -h for help
