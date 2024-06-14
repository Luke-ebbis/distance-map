#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Code adapted from
import Bio.PDB
import os
import itertools
import argparse
import numpy as np
from typing import List, Tuple, Set, Dict
from matplotlib import pyplot as plt

def load_structure(path: str):
    name = path
    if path.endswith("pdb"):
        structure = (Bio.PDB.PDBParser().get_structure(name, path))
    if path.endswith("cif"):
        structure = (Bio.PDB.MMCIFParser().get_structure(name, path))
    return structure


def calc_residue_dist(residue_one, residue_two):
    """Returns the C-alpha distance between two residues

    See [https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/
    protein_contact_map/]
    """
    # To get the side chains, you take the minimal distance between all the side
    # Chain atoms.
    try:
        diff_vector  = residue_one["CA"].coord - residue_two["CA"].coord
        return np.sqrt(np.sum(diff_vector * diff_vector))
    except KeyError as e:
        return np.nan


def calc_dist_matrix(chain_one, chain_two):
    """Returns a matrix of C-alpha distances between two chains

    See [https://warwick.ac.uk/fac/sci/moac/people/students/peter_cock/python/
    protein_contact_map/]
    """
    answer = np.zeros((len(chain_one), len(chain_two)), np.cfloat)
    for row, residue_one in enumerate(chain_one) :
        for col, residue_two in enumerate(chain_two) :
            answer[row, col] = calc_residue_dist(residue_one, residue_two)
    return answer


def contact_map(model,
                chain1: str,
                chain2: str,
                threshold=12.0):
    dist_matrix = calc_dist_matrix(model[chain1], model[chain2])
    contact_map = dist_matrix < threshold
    return contact_map

def generate_contact_maps(structure,
                          include: None | List[Set[str]] = None,
                          threshold=5.0,
                          folder=".",
                          names = Dict|None):
    os.makedirs(folder, exist_ok=True)
    print(f"Analysing distance maps at {threshold} Å")
    model = structure[0]
    if include:
        print(f"only investigating interactions between\n{include}")
    chains = [c.id for c in list(model.get_chains())]
    maps = dict()
    for pair in itertools.combinations(chains, 2):
        if include:
            pair = set(pair)
            if pair not in include:
                continue
        a, b = pair
        filename = f"contacts-{structure.id}-{a}-{b}-{threshold}.png"
        analyse_contacts(pair, structure, threshold,
                         f'{folder}/{filename}',
                         names)
    for pair in zip(chains, chains):
        if include:
            unique_pair = set(pair)
            if unique_pair not in include:
                continue
        if pair[0] == pair[1]:
            a, b = pair
            filename = f"contacts-{structure.id}-{a}-{b}-{threshold}.png"
            analyse_contacts((pair[0], pair[0]), structure, threshold,
                             f'{folder}/{filename}',
                             names)

def analyse_contacts(pair, structure, threshold, filename, names):
    a, b = sorted(list(pair))
    c = contact_map(structure[0],
                    a,
                    b,
                    threshold=threshold)
    # save only if sum is bigger than 1
    if c.sum() != 0:
        print(f"making contacts for {pair}\n {c.sum()} contacts in total")
        if names:
            a = names[a]
            b = names[b]
        make_contact_plot(c,
                          xname = a,
                          yname = b,
                          name=filename)

def make_contact_plot(contacts, xname, yname, name):
     fig = plt.figure()
     ax = fig.add_subplot(111)
     ax.imshow(np.transpose(contacts))
     ax.set_aspect('auto')
     plt.xlabel(f"Residue number on {xname}")
     plt.ylabel(f"Residue number on {yname}")
     fig.savefig(name)
     plt.close()

def parse_filter(x: str):
    if x:
        out = [set(item) for item in x.split(':')]
        print(f"filter: {out}")
        for item in out:
            print(len(item))
            if not len(item) <= 2:
                raise RuntimeError(f"You must give a pair of two symbols! Example 'ab' for chain a vs b.")
        return out

def parse_new_names(names: None|str):
    out = None
    if names:
        pairs = names.split(',')
        out = dict()
        for pair in pairs:
            key, value = pair.split(':')
            out[key] = value
    return out

def main():
    radius = 15
    parser = argparse.ArgumentParser(
                        prog='Contact map plotter',
                        description=f'Plotting contact maps {radius} Å.',
                        epilog='Inspired by the work of Peter Cock')
    parser.add_argument('structure',
                        help="A CIF or PDB file.")
    parser.add_argument('output',
                        help="Name of the directory to save the results")
    parser.add_argument('-f', '--filter',
                        default=None,
                        help="The chains to calculate specifically.\nGive each"
                        "pair as two symbols separated by :. For example"
                        "AB:CC")
    parser.add_argument('-n', "--names",
                        default=None,
                        help="A string of ''<id>:name' pairs separated by comma's.")
    parser.add_argument('-r', '--radius',
                        help="The radius for which two CA atoms are concidered in contact.")
    args = parser.parse_args()
    include = parse_filter(args.filter)
    new_names = parse_new_names(args.names)
    pdb_filename = args.structure
    structure = load_structure(pdb_filename)
    generate_contact_maps(structure,
                          include=include,
                          threshold=radius,
                          folder=args.output,
                          names=new_names)

if __name__ == "__main__":
    main()
