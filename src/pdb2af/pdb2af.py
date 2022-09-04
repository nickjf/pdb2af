#!/usr/bin/env python

import sys
import os
import argparse
import gzip
import json

import urllib.request

from importlib import resources
from pdbtools import pdb_selres

sys.tracebacklimit = 0
pdb2af_version = "0.0.20"

def renumber_atoms(pdb,out):
    pdb = open(pdb,'r')
    out = open(out,'w')
    
    count = 0
    for line in pdb:
        if 'ATOM' in line:
            count += 1
            out.write(line[0:6]+format(str(count)," >5s")+line[11:])
    out.write('TER\n')
    out.close()


def get_SIFTS(SIFTS_file="uniprot_segments_observed.tsv.gz"):
    data = urllib.request.urlretrieve('ftp://ftp.ebi.ac.uk/pub/databases/msd/sifts/flatfiles/tsv/'+SIFTS_file, SIFTS_file)

    pdb_chain_up_start_end = {}
    with gzip.open(SIFTS_file,'rb') as fin:        
        for line in fin:        
            line = line.decode("utf-8").split('\t')
            try:
                pdbid = line[0].upper()
                chain = line[1]
                upid = line[2]
                pdb_start = line[5]
                pdb_end = line[6]
                start = line[7]
                end = line[8].strip()

                if pdbid not in pdb_chain_up_start_end:
                    pdb_chain_up_start_end[pdbid] = {}
                if chain not in pdb_chain_up_start_end[pdbid]:
                    pdb_chain_up_start_end[pdbid][chain] = {}
                pdb_chain_up_start_end[pdbid][chain][upid] = {'pdb_start':pdb_start,'pdb_end':pdb_end,'start':start,'end':end}
            except:
                pass

    with resources.path("pdb2af.json","__init__.py") as f:
        SIFTS_path = os.path.dirname(f)+'/pdb_chain_up_start_end.json'
            
    with open(SIFTS_path, "w") as outfile:
        json.dump(pdb_chain_up_start_end, outfile)

    os.remove(SIFTS_file)


def main():

    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", type=str, help="file containing PDB IDs to match", required=True)
    parser.add_argument("-u", "--update", help="update SIFTS mapping between PDB & UniProt IDs", action="store_true")

    args = parser.parse_args()

    pdbs = []
    for line in open(args.input,'r'):
        if len(line.strip()) == 4:
            pdbs.append(line.strip().upper())

    if len(pdbs) == 0:
        print("-> ERROR did not find any PDB IDs in input file, quitting")
        sys.exit(1)
    else:

        if args.update:
            print("-> updating SIFTS file")
            get_SIFTS("pdb_chain_uniprot.tsv.gz")

        print("-> checking for SIFTS file")

        try:
            with resources.path("pdb2af.json","pdb_chain_up_start_end.json") as f:
                SIFTS_path = f
      
        except:

            print("-> couldn't find SIFTS file, attempting to download from PDBe")
            get_SIFTS("pdb_chain_uniprot.tsv.gz")


        try:
            with resources.path("pdb2af.json","pdb_chain_up_start_end.json") as f:
                SIFTS_path = f

            print("-> found SIFTS file")


            with open(SIFTS_path) as json_file:
                pdb_chain_up_start_end = json.load(json_file)

            for pdbid in pdbs:
                if pdbid in pdb_chain_up_start_end:
                    for chain in pdb_chain_up_start_end[pdbid]:
                        print("-> "+pdbid+" "+chain+" ",end='')
                        for upid in pdb_chain_up_start_end[pdbid][chain]:
                            print(" "+upid+" ",end='')
                            start = pdb_chain_up_start_end[pdbid][chain][upid]['start']
                            end = pdb_chain_up_start_end[pdbid][chain][upid]['end']
                            try:    
                                data = urllib.request.urlretrieve("https://alphafold.ebi.ac.uk/files/AF-"+upid+"-F1-model_v3.pdb", pdbid+"_"+chain+"_"+upid+".tmp")
                                _p,_r = pdb_selres.check_input(["-"+start+":"+end,pdbid+"_"+chain+"_"+upid+".tmp"])
                                new_pdb = pdb_selres.run(_p,sorted(_r))

                                out = open(pdbid+"_"+chain+"_"+upid+".tmp2",'w')
                                for line in enumerate(new_pdb):
                                    out.write(line[1])
                                
                                out.close()

                                renumber_atoms(pdbid+"_"+chain+"_"+upid+".tmp2", pdbid+"_"+chain+"_"+upid+".pdb")
                                
                                os.remove(pdbid+"_"+chain+"_"+upid+".tmp")
                                os.remove(pdbid+"_"+chain+"_"+upid+".tmp2")
                                print(" DONE")
                            except:
                                print(" NO MATCH")
                                pass

                else:
                    print("-> "+pdbid+" NO MATCH")


        except:
            print("-> ERORR couldn't find or download SIFTS file, quitting")
            sys.exit(1)





