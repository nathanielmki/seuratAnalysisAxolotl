import pandas as pd
import argparse
import sys, os
import csv
import numpy as np
import openpyxl

#def read_cluster_csv(csvfile):
    #print('read_cluster_csv(): type(csvfile)) = {}'.format(csvfile))
    #print('')
    #df = pd.read_csv(csvfile)
    #return df

def main():
    parser = argparse.ArgumentParser(description='Load csv for modification.')
    parser.add_argument('csvfile', type=argparse.FileType('r'), help='input csv file')
    args = parser.parse_args()

    #print('main(): type(args.csvfile)) = {}'.format(args.csvfile))
    #print('')

    df = pd.read_csv(args.csvfile)

    # look into groupby for pandas

    # df grouped by cluster column
    df1 = df.groupby(["cluster"])

    df2 = df1.apply(lambda x: x.sort_values(["avg_logFC"], ascending=False))
    #print(df2)

    # reset index
    df3 = df.reset_index(drop=True)
    df3.head()

    df_final = df3.groupby('cluster').head(1)
    print(df_final)

    df_final.to_excel("t1GenePerCluster_GeneIDs.xlsx")

if __name__ == '__main__':
    main()

