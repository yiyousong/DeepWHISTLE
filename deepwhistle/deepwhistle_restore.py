import argparse
import os
import sys
import tensorflow as tf
import numpy as np
import pandas as pd

tf.random.set_random_seed(7)
np.random.seed(7)
def bin2fasta(sequence):
    data_tmp = np.asarray(sequence)
    print('input data imported')
    data = []
    for i in range(data_tmp.shape[0]):
        tmp=''
        for j in range(1001):
            if data_tmp[i, j, 0] == 1:
                tmp+='A'
            elif data_tmp[i, j, 1] == 1:
                tmp += 'C'
            elif data_tmp[i, j, 2] == 1:
                tmp += 'G'
            elif data_tmp[i, j, 3] == 1:
                tmp += 'U'
            else:
                tmp += 'N'
        data.append(tmp)

    return data
def savefasta(data,loc):
 with open(loc,'w+') as w:
     for i in range(len(data)):
         w.write('\n>\n')
         w.write(data[i])
def read_fasta(input):
    fasta = []
    genename = []
    with open(input) as file_one:
        for line in file_one:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                active_sequence_name = line[1:]
                genename.append(active_sequence_name)
                continue
            sequence = line
            fasta.append(sequence)
    return genename, fasta


def process_long_sequence(genename, fasta):
    genename_1001 = []
    fasta_1001 = []
    for i in range(len(genename)):
        for j in range(len(fasta[i]) - 1000):
            genename_1001.append(genename[i] + ' position: ' + str(j))
            fasta_1001.append(fasta[i][j:j + 1001])
    return genename_1001, fasta_1001
def one_hot_encode(sequence):
    data_tmp = np.asarray(sequence)
    print('input data imported')
    data = np.zeros((data_tmp.shape[0], 1001, 4))
    for i in range(data_tmp.shape[0]):
        for j in range(1001):
            if data_tmp[i][j] == "A":
                data[i, j, :] = [1, 0, 0, 0]
            elif data_tmp[i][j] == "C":
                data[i, j, :] = [0, 1, 0, 0]
            elif data_tmp[i][j] == "G":
                data[i, j, :] = [0, 0, 1, 0]
            elif data_tmp[i][j] == "T":
                data[i, j, :] = [0, 0, 0, 1]
            elif data_tmp[i][j] == "U":
                data[i, j, :] = [0, 0, 0, 1]
            else:
                data[i, j, :] = [0, 0, 0, 0]
    return data


def restore(args):


    gen_input=np.asarray(pd.read_csv(args.gen_input))

    location = os.path.dirname(os.path.abspath(__file__))
    if not args.model_loc:
            if args.type == 'mature':
               model_loc='%s/model/model.h5'%(location)
            else:
                model_loc = '%s/model/full_model.h5' % (location)
    else:
        model_loc=args.model_loc
    model=tf.keras.models.load_model(model_loc)

    bp = 500


    if args.encoding == 'binary':
        data_tmp = pd.read_csv(args.input)
        data_tmp = np.asarray(data_tmp.fillna(0))
        data = np.reshape(data_tmp, (-1, 1001, 4))
        data=data[:,500-bp:501+bp,:]
        genename=[]

        prediction_final = model.predict([data,gen_input])
    else:
        genename_tmp, sequence_tmp = read_fasta(args.input)
        genename, sequence = process_long_sequence(genename_tmp, sequence_tmp)

        #check for DRACH motif
        sequence_predict = []
        genome_predict = []
        sequence_index = np.zeros(len(sequence))
        prediction_final = np.zeros(len(sequence))


        for i in range(len(sequence)):
            if sequence[i][500] == 'A':
                if (sequence[i][499] == 'G') | (sequence[i][499] == 'A'):
                    if (sequence[i][498] == 'A') | (sequence[i][498] == 'U') | (sequence[i][498] == 'T') | (
                            sequence[i][498] == 'G'):
                        if sequence[i][501] == 'C':
                            if (sequence[i][502] == 'A') | (sequence[i][502] == 'U') | (sequence[i][502] == 'T') | (
                                    sequence[i][502] == 'C'):
                                sequence_index[i] = 1
                                sequence_predict = np.append(sequence_predict, sequence[i])
                                genome_predict = np.append(genome_predict, gen_input[i])
        genome_predict=np.reshape(genome_predict,(-1,len(gen_input[0])))
        data=one_hot_encode(sequence_predict)
        data=data[:,500-bp:501+bp,:]

        prediction=model.predict([data,genome_predict])
        replace_index = 0
        for i in range(len(sequence_index)):
            if sequence_index[i] == 1:
                prediction_final[i] = prediction[replace_index]
                replace_index = replace_index + 1
        print('saving')
    return prediction_final,genename
def save_prediction(args):
    prediction_final, genename=restore(args)

    location = os.path.dirname(os.path.abspath(__file__))
    if args.encoding == 'binary':
        with open(args.output, "ab") as f:
            np.savetxt(f, prediction_final, fmt='%s',)
    else:
        if args.type=='mature':
            threshold = np.loadtxt('%s/m6A.txt' % (location))
        else:
            threshold = np.loadtxt('%s/m6A_full.txt' % (location))
        p_value = []

        for i in range(len(prediction_final)):
            p_value = np.append(p_value, sum(prediction_final[i] < threshold) / len(threshold))
        np.savetxt('%s/prediction.txt' % (args.output), prediction_final, fmt='%s')
        np.savetxt('%s/genename.txt' % (args.output), genename, fmt='%s')
        np.savetxt('%s/pvalue.txt' % (args.output), p_value, fmt='%s')
        print('finished')
if __name__ == '__main__':


    print('initializing preprocessing..')
    parser = argparse.ArgumentParser()
    parser.add_argument('-i','--input', default='fasta.txt', help='input file')
    parser.add_argument('-g','--gen_input', required=True, help='genome information input file')
    parser.add_argument('-o','--output', default='%s/output'%(os.path.dirname(os.path.abspath(__file__))), help='output folder')
    parser.add_argument('--model_loc', help='model location, please try default location first')
    parser.add_argument('--encoding', default='fasta', help='use default')
    parser.add_argument('-t','--type', default='mature', help='mature/full, default is mature')

    args = parser.parse_args()
    if args.encoding=='fasta':
        try:
            os.mkdir(args.output)
        except:
            print('unabe to create output folder, ignore this if output folder already exists')
        finally:
            save_prediction(args)
    else:
        save_prediction(args)
