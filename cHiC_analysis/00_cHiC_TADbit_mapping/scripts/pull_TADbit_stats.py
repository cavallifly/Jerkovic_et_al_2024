#!/usr/bin/env python
# coding: utf-8

import os
import sqlite3 as lite
from glob import glob
import pandas as pd
from regex import match

sample_dir = "/zdata/data/mdistefano/2021_04_13_Project_with_Ivana/cHiC_TADbit_mapping/"

samples = [f for f in os.listdir(sample_dir) if match('NPC', f)] + [f for f in os.listdir(sample_dir) if match('mESC', f)]
print(sample_dir + samples[0])
print(len(samples))

summary = {}
n = 0
multiple_names = {}
for sample in samples:
    n = n + 1
    sname = sample
    print(sample_dir, sample, "trace.db", n)
    con = lite.connect(os.path.join(sample_dir, sample, "trace.db"))
    summary[sname] = {
        'input total r1': 0, 
        'input total r2': 0, 
        'input RE-fragments r1': 0,
        'input RE-fragments r2': 0,
        'mapped full-reads r1': 0,
        'mapped full-reads r2': 0,
        'mapped RE-fragments r1': 0,
        'mapped RE-fragments r2': 0,
        'parsed reads r1 uniques': 0,
        'parsed reads r2 uniques': 0,
        'parsed reads r1 triple-wise': 0,
        'parsed reads r2 triple-wise': 0,
        'parsed reads r1 quadruple-wise': 0,
        'parsed reads r2 quadruple-wise': 0,
        'parsed reads r1 quintuple-wise': 0,
        'parsed reads r2 quintuple-wise': 0,
    }
    multiple_names[sname] = {
        "1": "triple-wise",
        "2": "quadruple-wise",
        "3": "quintuple-wise",
        "4": "sextuple-wise",
    }
    with con:
        cur = con.cursor()
        read_path = {}
        for read in [1, 2]:
            mids = {'full': [], 'frag': []}
            for mapping in ['full', 'frag']:
                cur.execute(("SELECT Entries, MAPPED_OUTPUTid FROM MAPPED_INPUTs "
                             "WHERE Read=%d and Frag='%s'") % (read, mapping))
                for m, mid in cur.fetchall():
                    if mapping=='full':
                        summary[sname][f'input total r{read}'] += m
                    elif mapping=='frag':
                        summary[sname][f'input RE-fragments r{read}'] += m
                    mids[mapping].append(mid)

            for mapping in ['full', 'frag']:
                for mid in mids[mapping]:
                    cur.execute('SELECT Uniquely_mapped, BEDid FROM MAPPED_OUTPUTs where PATHid=%s' % mid)
                    uniques, bedid = cur.fetchall()[0]  # get BEDid for parsed-outputs table
                    summary[sname][
                        f'mapped {"RE-fragments" if mapping == "frag" else "full-reads"} r{read}'
                    ] += uniques
                    read_path[read] = bedid

            cur.execute('SELECT Total_uniquely_mapped, Multiples FROM PARSED_OUTPUTs where PATHid=%s' %
                        read_path[read])

            uniques, multiples = cur.fetchall()[0]
            multiples = dict((v.split(':')[0], int(v.split(':')[1])) for v in multiples.split(','))
            summary[sname][f'parsed reads r{read} uniques'] = uniques
            for k in multiples:
                summary[sname][f'parsed reads r{read} {multiple_names[sname][k]}'] = multiples[k]


            tot = sum(summary[sname][k] for k in summary[sname] 
                      if k.startswith('mapped') and k.endswith(f'r{read}'))
            cur.execute('SELECT Multiples FROM PARSED_OUTPUTs where Total_uniquely_mapped=%s' % tot)
            m = cur.fetchall()[0][0]

        cur.execute('SELECT Total_interactions, Multiple_interactions FROM INTERSECTION_OUTPUTs')
        count, mult = cur.fetchall()[0]
        for elt in mult.split():
            m, c = map(int, elt.split(':'))
            count += (m + 1) * m / 2 * c - c
        summary[sname]['intersection'] = count
        cur.execute('SELECT Name, Count, Applied FROM FILTER_OUTPUTs')
        for name, num, applied in cur.fetchall():
            if applied == "True":
                summary[sname][f"Filtered {name}"] = num
            elif applied == "False":
                summary[sname][f"Un-filtered {name}"] = num
            else:
                summary[sname][f"{name}"] = num

df = pd.DataFrame.from_dict(summary)
df.to_excel('HiC-processing_stats.xls', merge_cells=False, 
            sheet_name="processing_stats")

out = open('TADbit_mapping_stats_per_sample.txt', 'w')
out.write(f"Sample\tInput_read_pairs\tMapped_full_r1\tMapped_frag_r1\tMapped_full_r2\tUniquely_mapped_pairs\tValid-pairs\n")
for item in summary.items():
    out.write("%-20s\t%-15d\t%-15d\t%-15d\t%-15d\t%-15d\t%-15d\t%-15d\n" %(item[0],item[1]['input total r1'],item[1]['mapped full-reads r1'],item[1]['mapped RE-fragments r1'],item[1]['mapped full-reads r2'],item[1]['mapped RE-fragments r2'],item[1]['intersection'],item[1]['valid-pairs']))
out.close
