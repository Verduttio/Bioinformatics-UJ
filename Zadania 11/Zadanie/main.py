import numpy as np
import pandas as pd
import Bio.Cluster


def check_data_presence(value):
    return 0 if value == 'x' else 1


def convert_or_zero(value):
    try:
        return float(value)
    except ValueError:
        return 0.0


if __name__ == '__main__':
    file_name = 'yeast_expression.txt'
    yeast_data_raw = pd.read_csv(file_name, sep='\t', header=None, skiprows=1,
                                 usecols=[0, 1, 2, 3, 4, 5, 6]).transpose()
    yest_mask = np.array([[check_data_presence(x) for x in row] for row in yeast_data_raw.values])
    yeast_data_converted = np.array([[convert_or_zero(x) for x in row] for row in yeast_data_raw.values])

    clusterid, error, nfound = Bio.Cluster.kcluster(yeast_data_converted, nclusters=2, mask=yest_mask, transpose=False,
                                                    npass=50, method='a', dist='c')

    print(clusterid)
