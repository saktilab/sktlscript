#!/usr/bin/env python3

import sys
import argparse

parser = argparse.ArgumentParser(description="Program untuk menjalankan CREST")
parser.add_argument('-m','--method',type=str, default="gfnff", help="Metode yang digunakan untuk menjalankan program CREST. Default: gfnff. Pilihan: gfn1, gfn2")
parser.add_argument('-i','--input',type=str, default="struc.xyz",help="Input struktur yang digunakan. Default: struc.xyz.")
parser.add_argument('-xtbin','--xtbin',type=str,default='/home/coto/software/xtb/build/xtb', help="Lokasi binary xtb yang digunakan.")
parser.add_argument('-crestbin','--crestbin',type=str,default='/home/adit/opt/crest/crest')
parser.add_argument('-mode','--mode',type=str, default='none', help="Jenis perhitungan yang dilakukan. Tidak ada nilai default, pilihan: nci, proton, tautomer, cregen, mdopt")
parser.add_argument('-ewin','--ewin', type=float, default=2.0, help="Besarnya jendela energi yang digunakan untuk pencarian konformasi. Default: 2.0 kcal/mol.")
parser.add_argument('-gbsa','--gbsa', type=str, help="Model pelarut GBSA yang digunakan.")
parser.add_argument('-alpb','--alpb', type=str, help="Model pelarut ALPB yang digunakan.")
parser.add_argument('-np','--nproc',type=int, default=1, help='Jumlah processor yang digunakan untuk menghitung.')
parser.add_argument('-cregen','--cregen',type=str, default='xtb.trj', help='File trayektori yang digunakan untuk seleksi struktur manual menggunakan CREST.')
parser.add_argument('-mdopt','--mdopt',type=str, default='xtb.trj', help='File trayektori dinamika molekul yang ingin dioptimasi secara manual.')

opt = parser.parse_args(sys.argv[1:])

if opt.mode == 'none':
	os.system('{} {} --{} -ewin {} -xnam {} ')
elif opt.mode == 'cregen':

elif opt.mode == 'mdopt':

else:
