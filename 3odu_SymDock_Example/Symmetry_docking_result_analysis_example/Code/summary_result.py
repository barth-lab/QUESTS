#! /usr/bin/env python
import sys, shutil

def main():
    filename = sys.argv[1]
    fileinfor = readin_cross_angle(filename)
#    filename_list = [i[1] for i in fileinfor]
#    for filename in filename_list:
#        shutil.copy(filename+".pdb", "../selected_structure/")

def readin_cross_angle(filename):
    infor = []
    lineinfor = file(filename, 'r').readlines()
    for line in lineinfor:
        Isc, name, a, cross_angle, b, distance = line.strip().split()
        if float(cross_angle) <= 85 and float(cross_angle) >= -50 and float(distance) <= 40:
#            infor.append((name, float(cross_angle), float(distance)))
             print line.strip()
        else: continue
    return infor

if __name__ == '__main__':
    main()
