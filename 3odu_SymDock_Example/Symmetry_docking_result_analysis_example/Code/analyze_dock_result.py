#! /usr/bin/env python

import numpy as np
from Structure_tools import *
from Vector_tools import *
import sys,math
#from chempy import cpv
import cpv 

def main():
    dock_result = sys.argv[1]
    Isc, dock_result = dock_result.split(",")
    resi, resj = sys.argv[2], sys.argv[3]
    resx, resy = sys.argv[4], sys.argv[5]
    prostr = Protein_Structure(dock_result )
    v1, m1, raw_data1 = fitting_helical_axis(resi, resj, prostr)
    v2, m2, raw_data2 = fitting_helical_axis(resx, resy, prostr)
    #print "v1,v2", tuple(v1), m1, tuple(v2), m2

#    a1 = math.fsum(t(v1*v2))
#    a2 = vl1*vl2
#    final_angle = math.acos( a1/a2 )
    final_angle = Np_Angle_Between_Vector(v1,v2, length(v1), length(v2))
    # please note that this is not the closest point between two helices.
    # This is just for compare cross angle the of sampe pair
    # I am just measure the middle point of the same helices,
    final_distance = DistanceMeasure(m1, m2)
    v3 = np.cross(v1, v2)
    AB_vector = Vector_generate(m1, m2)
    AB_vector = np.array(AB_vector)
    #print "v3, ab_vecotr", v3, AB_vector, vl3, final_distance
    direction_angle = Np_Angle_Between_Vector(v3, AB_vector, length(v3), final_distance)
    if direction_angle <= math.pi/2 : final_angle = -final_angle
    #print direction_angle/math.pi*180, final_angle/math.pi*180, final_distance
    print "%s\t%s\t cross-angle  %.2f\tdistance  %.2f" %(Isc, dock_result, final_angle[0]/math.pi*180, final_distance)


    # plot the line
#    line1 = v1 * np.mgrid[-10:10:2j][:, np.newaxis]
#    line1 += m1
#    line2 = v2 * np.mgrid[-10:10:2j][:, np.newaxis]
#    line2 += m2

#    import matplotlib.pyplot as plt
#    import mpl_toolkits.mplot3d as m3d

#    ax = m3d.Axes3D(plt.figure())
#    ax.scatter3D(*raw_data1.T)
#    ax.scatter3D(*raw_data2.T)
#    ax.plot3D(*line1.T)
#    ax.plot3D(*line2.T)
#    plt.show()


def fitting_helical_axis(resa, resb, prostr):
    helical_points = Readin_helical_region(resa, resb, prostr)
    data_mean = helical_points.mean(axis=0)
    uu, s, vh = np.linalg.svd(helical_points - data_mean)
    vec = cpv.normalize(vh[0])
    if dotproduct(vec, helical_points[-1] - helical_points[0]) < 0:
        vec = cpv.negate(vec)
    # vv[0] is the direction
    # v = vv.conj().transpose()
    return vec, data_mean, helical_points


def Readin_helical_region(resi, resj, prostr):
    collection = []
    for res in sorted(prostr.Seq_coord.keys()):
        if res[0] >= int(resi) and res[0] <= int(resj):
            for Atom in prostr.Seq_coord[res]:
                    if Atom.AtomName in [ "CA" ]:
                        collection.append( Atom.AtomLocation )
                    else: continue
        else: continue
    collection = np.array(collection)
    return collection



if __name__ == '__main__':
    main()
