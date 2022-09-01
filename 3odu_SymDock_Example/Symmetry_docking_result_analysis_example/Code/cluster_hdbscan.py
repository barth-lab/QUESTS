#! /usr/bin/env python
# used to cluster the conformations of the CXCR4 dimer
# general used parameters:  

import sys, hdbscan
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sns.set_context('poster')
sns.set_color_codes()
plot_kwds = {'alpha' : 0.25, 's' : 80, 'linewidths':0}


def plot_clusters(data, clusterer):
    labels = clusterer.fit_predict(data)
    palette = sns.color_palette('deep', np.unique(labels).max() + 1)
    colors = [palette[x] if x >= 0 else (0.0, 0.0, 0.0) for x in labels]
    plt.scatter(data.T[0], data.T[1], c=colors, **plot_kwds)
    frame = plt.gca()
    frame.axes.get_xaxis().set_visible(True)
    frame.axes.get_yaxis().set_visible(True)
    print "showing plot"
    plt.show()


def main():
    "test a version a scale factor of 5.0 is considered for the angle difference"
    # print sys.argv
    datafile = sys.argv[1]
    # energy filter index d: when d = 10 indicate that we only analyze the top 10 percentage of the
    # data set by interface energy
    d = int(sys.argv[2])
    data_array = data_to_array(datafile)
    new_array = data_array[Energy_filter(data_array,d)]
    #new_array = data_array
    #print new_array[3,:]
    data_infor = new_array[:, 1:]
    angle_list = new_array[:, 1]
    distance_list = 5*new_array[:, 2]
    # reassign the distance value with a factor of 5 to compensate for small sampling space
    for i in range(len(angle_list)):
        data_infor[i,0] = angle_list[i]
        data_infor[i,1] = distance_list[i]
    # print distance_list
    # print len(data_infor)
    new_array = data_array[Energy_filter(data_array,d)]
    clusterer = hdbscan.HDBSCAN(min_cluster_size=4, metric="cityblock")
    # print clusterer.metric
    clusterer.fit_predict(data_infor)
    #print cluster_labels
    #print len(clusterer.labels_)
    # check number of the clusters
    m = clusterer.labels_.max()
    #print m
    for i in range(0,m+1):
        print len(data_infor[clusterer.labels_==i]), new_array[clusterer.labels_==i][0]
        # print new_array[clusterer.labels_==i]
        output_clustering_result(new_array[clusterer.labels_==i],i)
    plot_clusters(data_infor, clusterer)


def data_to_array(datafile):
    data = np.fromfile(datafile, sep=" ")
    n = len(data)
    data_array = np.resize(data, (n/3, 3))
    #print data_array
    return data_array


def Energy_filter(data_array, d):
    energy_list = data_array[:, 0]
    n = len(energy_list)
    if d > 1:
      filter = sorted(energy_list)[n/d]
      print sorted(energy_list)[n/d]
    else:
        filter = sorted(energy_list)[n-1]
        print sorted(energy_list)[n-1]
    return energy_list <= filter


def output_clustering_result(new_array, i):
    outputfile = file('cluster%d.txt' %i,'w')
    A,B = new_array.shape
    print A, B
    for a in range(A):
        print >> outputfile, "%s\t%s\t%s" %(tuple(new_array[a,:]))
    print >> outputfile, "Average: %.2f\tSize: %d" %(sum(new_array[:,0])/A, A)
    angle_min, angle_max = sorted(new_array[:,1])[0], sorted(new_array[:,1])[-1]
    dist_min, dist_max = sorted(new_array[:,2])[0], sorted(new_array[:,2])[-1]
    print >> outputfile, "angle_min: %.2f\tangle_max: %.2f" %(angle_min, angle_max)
    print >> outputfile, "dist_min: %.2f\tdist_max: %.2f" %(dist_min, dist_max)

if __name__=="__main__":
    main()
