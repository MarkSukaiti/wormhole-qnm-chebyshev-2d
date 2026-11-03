#!/usr/bin/env python3
"""Reference values for the Kerr benchmark (Sec. V of the paper).

Massless scalar (s = 0) Kerr quasinormal frequencies M*omega from Leaver's continued-fraction
method as implemented in the qnm package (L. C. Stein, J. Open Source Softw. 4, 1683 (2019);
https://github.com/duetosymmetry/qnm), for l = 1..9, m = +1 and -1, overtones n = 0..11, at
the spins a/M = 0, 0.1, ..., 0.9, 0.99 of the benchmark.  qnm follows every overtone
continuously in the spin from the Schwarzschild limit, and its labels are adopted.

This is the only Python script of the repository: the qnm package is a Python package.  Its
output, kerr_qnm_reference.csv, is included, so it only needs to be run to regenerate it.

Usage:   pip install qnm ;  python3 qnm_reference.py      (about an hour on two cores)
"""
import csv
import multiprocessing as mp
import warnings

import qnm

warnings.filterwarnings('ignore')
SPINS = [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
KEYS = [(l, m, n) for l in range(1, 10) for m in (1, -1) for n in range(12)]


def job(key):
    l, m, n = key
    try:
        seq = qnm.modes_cache(s=0, l=l, m=m, n=n)
        return key, [complex(seq(a=a)[0]) for a in SPINS]
    except Exception as err:             # a sequence that qnm cannot follow is skipped
        return key, str(err)


if __name__ == '__main__':
    qnm.download_data()                  # precomputed Schwarzschild seeds
    with mp.Pool() as pool:
        res = dict(pool.imap_unordered(job, KEYS))
    with open('kerr_qnm_reference.csv', 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['l', 'm', 'n', 'a_over_M', 'Re_Momega', 'Im_Momega'])
        for key in sorted(res):
            if isinstance(res[key], str):
                print('skipped', key, res[key])
                continue
            for a, om in zip(SPINS, res[key]):
                w.writerow([*key, a, repr(om.real), repr(om.imag)])
    print('written kerr_qnm_reference.csv')
