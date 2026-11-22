"""Step C1: our HKY ts/tv parameterisation vs p4's 2p matrix.

The kappa behind a given ts/tv ratio is not a constant, so both sides are driven
from the same pi and the same tstv. The quantity compared is the log-likelihood
of one fixed four-tip alignment, so the reported difference is a scalar per
(pi, kappa) pair, not an entry-wise difference of the Q matrices themselves.
"""
import io
import json
import os
import subprocess
import sys

import numpy as np

from _env import directory

W = directory('SEQ_SIM_EXT_WORK', 'a writable working directory')
P4 = directory('P4_DIR', 'the built p4-phylogenetics checkout')
SRC = os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))), "src")
sys.path.insert(0, SRC)
from seq_sim.core.models import setup_model
from seq_sim.tree.tree_parser import parse_trees

B = os.path.join(W, 'b2')
os.makedirs(B, exist_ok=True)
TREE = '((A:0.1,B:0.1):0.2,(C:0.1,D:0.1):0.2);'
pats = [('A','A','A','A'), ('A','G','A','G'), ('G','G','G','G'), ('A','C','G','T')]
REPS = 250
rows = {'A': [], 'B': [], 'C': [], 'D': []}
for p in pats:
    for _ in range(REPS):
        for k in range(4):
            rows['ABCD'[k]].append(p[k])
seqs = {k: ''.join(v) for k, v in rows.items()}
n = len(seqs['A'])
nex = ['#NEXUS', 'BEGIN DATA;']
nex.append('  DIMENSIONS NTAX=4 NCHAR=%d;' % n)
nex += ['  FORMAT DATATYPE=DNA MISSING=? GAP=-;', '  MATRIX']
for k in 'ABCD':
    nex.append('  %s %s' % (k, seqs[k]))
nex += ['  ;', 'END;', 'BEGIN TREES;', '  TREE t1 = ' + TREE, 'END;']
io.open(os.path.join(B, 'dna.nex'), 'w', encoding='utf-8').write("\n".join(nex) + "\n")
cfgs = [{'tag': 'equal_k8', 'pi': [0.25,0.25,0.25,0.25], 'kappa': 8.0},
        {'tag': 'skew_k6p25', 'pi': [0.4,0.1,0.1,0.4], 'kappa': 6.25},
        {'tag': 'mid_k5', 'pi': [0.3,0.2,0.2,0.3], 'kappa': 5.0}]
io.open(os.path.join(B, 'hky.json'), 'w').write(json.dumps(cfgs))
m0 = setup_model('HKY', tstv=2.0, freqs=[0.25]*4)
print('state_chars', m0.state_chars)
def ours(pi, kappa):
    r = pi[0] + pi[2]; y = pi[1] + pi[3]; h = pi[0]*pi[2] + pi[1]*pi[3]
    tstv = kappa * h / (r * y)
    m = setup_model('HKY', tstv=tstv, freqs=pi, equal_exchangeability=False)
    t = parse_trees(TREE)[0][2]
    st = {c: i for i, c in enumerate(m.state_chars)}
    def Pmat(ch):
        buf = np.zeros(16); m.set_matrix(buf, ch.length0 * ch.branch_rate)
        M = buf.reshape(4,4); P = np.empty((4,4)); P[:,0] = M[:,0]; P[:,1:] = np.diff(M, axis=1); return P
    tot = 0.0
    for p in pats:
        obs = {}
        for j, c in enumerate('ABCD'):
            obs[c] = p[j]
        def nl2(nd):
            if nd.tip_no >= 0:
                v = np.zeros(4); v[st[obs[nd.name]]] = 1.0; return v
            L = np.ones(4)
            for ch in (nd.branch1, nd.branch2):
                L = L * (Pmat(ch) @ nl2(ch))
            return L
        tot += REPS * np.log(float(np.array(m.freq) @ nl2(t.root)))
    return tot, tstv
print('%-12s %22s %22s %10s %10s' % ('tag', 'p4', 'ours', 'diff', 'tstv'))
for i, c in enumerate(cfgs):
    out = subprocess.run(['python3', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'p4_hky_side.py'), str(i)], cwd=P4, capture_output=True, text=True, env=dict(os.environ, PYTHONPATH=P4))
    line = [l for l in out.stdout.splitlines() if l.startswith('RES')]
    if not line:
        raise SystemExit('p4 produced no RES line')
    p4v = float(line[0].split()[2])
    ov, tstv = ours(c['pi'], c['kappa'])
    print('%-12s %22.12f %22.12f %10.1e %10.6f' % (c['tag'], p4v, ov, p4v - ov, tstv))
