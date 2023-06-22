#!/usr/bin/env python
#-*- coding: utf-8 -*-

from argparse import ArgumentParser
from pathlib import Path

### Get arguments ###
def get_option():
    topf = "topol.top"
    psff = "out.psf"
    argparser = ArgumentParser()
    argparser.add_argument('-f', metavar='top_file', type=str,
                            default=topf,
                            help='input top file name (default: topol.top)')
    argparser.add_argument('-o', metavar='psf_file', type=str,
                            default=psff,
                            help='output psf file name (default: out.psf)')
    return argparser.parse_args()


### Read a top file and check mols and itps ###
def read_top(fname):
    mol_lst = []
    itps = []
    flag = None
    dir_path = Path(fname).parents[0]
    with open(fname, "r") as f:
        for line in f.readlines():
            line = line.strip()
            data = line.split()
            if len(data) == 0 or line.startswith(';'):
                continue
            if line.startswith('#include'):
                itps.append(dir_path / data[1].replace('"', ''))
                continue
            if line.startswith('['):
                flag = data[1]
                continue
            if flag == "molecules":
                resname = data[0]
                mol_lst.append([resname, int(data[1])])
    return mol_lst, itps


### Output the top file info ###
def output_info(mol_lst, itps):
    print("System information")
    for key, val in mol_lst:
        print("-", key, "#", val)
    print("Included itps")
    for itp in itps:
        print("-", itp)


### Read included itps and make molecule/atom dict. ###
def read_itp(itps):
    mol_dct = {}
    atm_dct = {}
    for itp in itps:
        try:
            f = open(itp, "r")
        except IOError:
            print(f"WARNING: {itp} cannot be opened") 
            continue
        for line in f.readlines():
            line = line.strip()
            iskip = line.find(";")
            if iskip != -1:
                line = line[:iskip]
            data = line.split()
            if len(data) == 0 or line.startswith('#'):
                continue
            if line.startswith('['):
                read_flag = data[1]
                continue
            if read_flag == 'moleculetype':
                molname = data[0]
                mol_dct[molname] = {}
                mol_dct[molname]["atom"] = [] # [atomtype, resname, atomname, charge, mass, massflag]
                mol_dct[molname]["bond"] = [] # [b1, b2]
                mol_dct[molname]["angle"] = [] # [a1, a2, a3]
                mol_dct[molname]["dihed"] = [] # [d1, d2, d3, d4]
                mol_dct[molname]["impro"] = [] # [i1, i2, i3, i4]
            elif read_flag == 'atoms':
                atype = data[1]
                rname = data[3]
                aname = data[4]
                charge = data[6]
                if len(data) == 8:
                    mass = data[7]
                    mflag = None
                else:
                    mass = "10.0"
                    mflag = True
                mol_dct[molname]["atom"].append([atype, rname, aname, charge, mass, mflag])
            elif read_flag == 'bonds' or read_flag == 'constraints':
                b1 = data[0]
                b2 = data[1]
                mol_dct[molname]["bond"].append([b1, b2])
            elif read_flag == 'angles':
                a1 = data[0]
                a2 = data[1]
                a3 = data[2]
                mol_dct[molname]["angle"].append([a1, a2, a3])
            elif read_flag == 'settles':
                mol_dct[molname]["bond"].append(['1', '2'])
                mol_dct[molname]["bond"].append(['1', '3'])
                #mol_dct[molname]["bond"].append(['2', '3'])
                mol_dct[molname]["angle"].append(['2', '1', '3'])
            elif read_flag == 'dihedrals':
                dihed_type = data[4]
                if dihed_type in ['1', '3', '5', '8', '9']:
                    d1 = data[0]
                    d2 = data[1]
                    d3 = data[2]
                    d4 = data[3]
                    mol_dct[molname]["dihed"].append([d1, d2, d3, d4])
                elif dihed_type in ['2', '4']:
                    i1 = data[0]
                    i2 = data[1]
                    i3 = data[2]
                    i4 = data[3]
                    mol_dct[molname]["impro"].append([i1, i2, i3, i4])
                else:
                    print(f"WARNING: Dihedral func. type {dihed_type} cannot be assigned") 
            elif read_flag == 'atomtypes':
                atype = data[0]
                if len(data) < 7:
                    mass = data[1]
                    charge = data[2]
                else:
                    mass = data[2]
                    charge = data[3]
                atm_dct[atype] = [mass, charge]
            else:
                continue
        f.close()
    return mol_dct, atm_dct


### Write psf ##
def write_psf(fname, mol_lst, mol_dct, atm_dct):
    ncomments = 2
    ndon = nacc = nnb = 0
    natoms = nbonds = nangles = ndiheds = nimpros = 0
    for mol, num in mol_lst:
        natoms += int(len(mol_dct[mol]["atom"])) * num
        nbonds += int(len(mol_dct[mol]["bond"])) * num
        nangles += int(len(mol_dct[mol]["angle"])) * num
        ndiheds += int(len(mol_dct[mol]["dihed"])) * num
        nimpros += int(len(mol_dct[mol]["impro"])) * num
    fpsf = open(fname, "w")
    print('PSF', file=fpsf)
    print(file=fpsf)
    print("{:8} !NTITLE".format(ncomments), file=fpsf)
    print("*", file=fpsf)
    print("*", file=fpsf)
    print(file=fpsf)
    ## Atom
    print("{:8} !NATOM".format(natoms), file=fpsf)
    unuse = 0
    nid = 0
    wflag = True
    for mol, num in mol_lst:
        segid = 0
        for ir in range(num):
            segid += 1
            for ia in range(len(mol_dct[mol]["atom"])):
                nid += 1
                atype = mol_dct[mol]["atom"][ia][0]
                rname = mol_dct[mol]["atom"][ia][1]
                aname = mol_dct[mol]["atom"][ia][2]
                charge = mol_dct[mol]["atom"][ia][3]
                mass = mol_dct[mol]["atom"][ia][4]
                mflag = mol_dct[mol]["atom"][ia][5]
                if len(mol) > 4:
                    mname = mol[:4]
                else:
                    mname = mol
                if mflag:
                    try:
                       mass = atm_dct[atype][0]
                    except KeyError:
                       if wflag:
                           print("WARNING: Mass values are not correctly assinged from itps")
                           wflag = False
                fpsf.write("%8d %-5s%-8d%-5s%-5s%-5s%12.6f%12.4f%12d\n" %
                           (nid, mname, segid, rname, aname, atype,
                            float(charge), float(mass), unuse))
    print(file=fpsf)
    ## Bond
    print("{:8} !NBOND".format(nbonds), file=fpsf)
    cnt = 0
    atmtot = 0
    for mol, num in mol_lst:
        for ir in range(num):
            for ib in range(len(mol_dct[mol]["bond"])):
                fpsf.write("%8d%8d" % (int(mol_dct[mol]["bond"][ib][0]) + atmtot,
                                       int(mol_dct[mol]["bond"][ib][1]) + atmtot))
                cnt += 1
                if cnt % 4 == 0:
                    fpsf.write("\n")
            atmtot += len(mol_dct[mol]["atom"])
    if cnt % 4 != 0:
        fpsf.write("\n")
    print(file=fpsf)
    ## Angle
    print("{:8} !NTHETA".format(nangles), file=fpsf)
    cnt = 0
    atmtot = 0
    for mol, num in mol_lst:
        for ir in range(num):
            for ia in range(len(mol_dct[mol]["angle"])):
                fpsf.write("%8d%8d%8d" % (int(mol_dct[mol]["angle"][ia][0]) + atmtot,
                                          int(mol_dct[mol]["angle"][ia][1]) + atmtot,
                                          int(mol_dct[mol]["angle"][ia][2]) + atmtot))
                cnt += 1
                if cnt % 3 == 0:
                    fpsf.write("\n")
            atmtot += len(mol_dct[mol]["atom"])
    if cnt % 3 != 0:
        fpsf.write("\n")
    print(file=fpsf)
    ## Proper dihedral
    print("{:8} !NPHI".format(ndiheds), file=fpsf)
    cnt = 0
    atmtot = 0
    for mol, num in mol_lst:
        for ir in range(num):
            for ip in range(len(mol_dct[mol]["dihed"])):
                fpsf.write("%8d%8d%8d%8d" %
                           (int(mol_dct[mol]["dihed"][ip][0]) + atmtot,
                            int(mol_dct[mol]["dihed"][ip][1]) + atmtot,
                            int(mol_dct[mol]["dihed"][ip][2]) + atmtot,
                            int(mol_dct[mol]["dihed"][ip][3]) + atmtot))
                cnt += 1
                if cnt % 2 == 0:
                    fpsf.write("\n")
            atmtot += len(mol_dct[mol]["atom"])
    if cnt % 2 != 0:
        fpsf.write("\n")
    print(file=fpsf)
    ## Improper dihedral
    print("{:8} !NIMPHI".format(nimpros), file=fpsf)
    cnt = 0
    atmtot = 0
    for mol, num in mol_lst:
        for ir in range(num):
            for ii in range(len(mol_dct[mol]["impro"])):
                fpsf.write("%8d%8d%8d%8d" %
                           (int(mol_dct[mol]["impro"][ii][0]) + atmtot,
                            int(mol_dct[mol]["impro"][ii][1]) + atmtot,
                            int(mol_dct[mol]["impro"][ii][2]) + atmtot,
                            int(mol_dct[mol]["impro"][ii][3]) + atmtot))
                cnt += 1
                if cnt % 2 == 0:
                    fpsf.write("\n")
            atmtot += len(mol_dct[mol]["atom"])
    if cnt % 2 != 0:
        fpsf.write("\n")
    print(file=fpsf)
    ## Others
    print("{:8} !NDON".format(ndon), file=fpsf)
    print(file=fpsf)
    print("{:8} !NACC".format(nacc), file=fpsf)
    print(file=fpsf)
    print("{:8} !NNB".format(nnb), file=fpsf)
    print(file=fpsf)
    print()
    print("Done")
    fpsf.close()


def main():
    print("top2psf")
    args = get_option()
    topf = args.f
    psff = args.o
    mol_lst, itps = read_top(topf)
    output_info(mol_lst, itps)
    mol_dct, atm_dct = read_itp(itps)
    write_psf(psff, mol_lst, mol_dct, atm_dct)


if __name__ == '__main__':
    main()
