"""Step C2: p4 setModelComponentOnNode(clade=1) vs our per-edge composition.

The question is which edges p4 actually transforms when a clade is named, and
whether our [f=N] annotation covers the same edge set.
"""
import io
import itertools
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
from seq_sim.core.aa_models import AA_ORDER, AAModel, empirical_rates
from seq_sim.core.gamma import discrete_gamma
from seq_sim.core.profile_mixture import branch_transform, build_profile_models
from seq_sim.tree.tree_parser import parse_trees

B = os.path.join(W, 'b2')
os.makedirs(B, exist_ok=True)
TREE = '((A:0.1,B:0.1)H1:0.2,(C:0.1,D:0.1)H2:0.2);'
pats = [('A','A','A','A'), ('A','G','A','G'), ('G','G','G','G'), ('A','R','N','D')]
REPS = 250
K = 4
ALPHA = 0.5
GW, CR = discrete_gamma(ALPHA, K)
print('cat rates', np.round(CR, 8))
rows = {'A': [], 'B': [], 'C': [], 'D': []}
for p in pats:
    for _ in range(REPS):
        for k in range(4):
            rows['ABCD'[k]].append(p[k])
seqs = {k: ''.join(v) for k, v in rows.items()}
n = len(seqs['A'])
nex = ['#NEXUS', 'BEGIN DATA;']
nex.append('  DIMENSIONS NTAX=4 NCHAR=%d;' % n)
nex += ['  FORMAT DATATYPE=PROTEIN MISSING=? GAP=-;', '  MATRIX']
for k in 'ABCD':
    nex.append('  %s %s' % (k, seqs[k]))
nex += ['  ;', 'END;', 'BEGIN TREES;', '  TREE t1 = ' + TREE, 'END;']
io.open(os.path.join(B, 'aa.nex'), 'w', encoding='utf-8').write("\n".join(nex) + "\n")
base, _pn, _wt = build_profile_models('C10', empirical_rates('LG'), {})
pi0 = np.array(base[0].freq)
pit = np.array(branch_transform(pi0, 3.0, 0.4))
io.open(os.path.join(B, 'clade.json'), 'w').write(json.dumps([{'pi0': list(pi0), 'pit': list(pit)}]))
m0 = AAModel(pi0, empirical_rates('LG'))
m1 = AAModel(pit, empirical_rates('LG'))
st = {c: i for i, c in enumerate(AA_ORDER)}
def ours(targets):
    t = parse_trees(TREE)[0][2]
    P = {}
    for ch in t.node_list:
        nm = ch.name
        m = m1 if nm in targets else m0
        for k in range(K):
            buf = np.zeros(400)
            pr = ch.branch_rate[0] if isinstance(ch.branch_rate, (list, tuple)) else ch.branch_rate
            m.set_matrix(buf, ch.length0 * pr * CR[k])
            M = buf.reshape(20,20)
            Q = np.empty((20,20)); Q[:,0] = M[:,0]; Q[:,1:] = np.diff(M, axis=1)
            P[(id(ch), k)] = Q
    def nl(nd, obs, k):
        if nd.tip_no >= 0:
            v = np.zeros(20); v[st[obs[nd.name]]] = 1.0; return v
        L = np.ones(20)
        for ch in (nd.branch1, nd.branch2):
            L = L * (P[(id(ch), k)] @ nl(ch, obs, k))
        return L
    prior = pit if 'root' in targets else pi0
    tot = 0.0
    for p in pats:
        obs = dict(zip('ABCD', p))
        acc = 0.0
        for k in range(K):
            acc += GW[k] * float(prior @ nl(t.root, obs, k))
        tot += REPS * np.log(acc)
    return tot
ANCHOR = -7573.131275040946
a = ours({'A'})
print('anchor ours(A)  : %.12f   diff %+.2e' % (a, a - ANCHOR))
out = subprocess.run(['python3', os.path.join(os.path.dirname(os.path.abspath(__file__)), 'p4_clade_side.py')], cwd=P4, capture_output=True, text=True, env=dict(os.environ, PYTHONPATH=P4))
line = [l for l in out.stdout.splitlines() if l.startswith('RES')]
if not line:
    print('P4 FAIL', out.stderr[-400:]); sys.exit(1)
p4v = float(line[0].split()[1])
print('p4 clade=1 on H1 : %.12f' % p4v)
print('ours {H1}        : %.12f  diff %+.3e' % (ours({'H1'}), ours({'H1'}) - p4v))
slots = ['root','H1','H2','A','B','C','D']
hits = []
for r in range(8):
    for comb in itertools.combinations(slots, r):
        v = ours(set(comb))
        if abs(v - p4v) < 1e-4:
            hits.append((comb, v))
print('HITS', len(hits), hits[:4])
