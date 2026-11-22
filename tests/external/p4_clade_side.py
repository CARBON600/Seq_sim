import json
import os

import p4

from _env import directory

W = directory('SEQ_SIM_EXT_WORK', 'a writable working directory')
B = os.path.join(W, 'b2')
cfg = json.load(open(os.path.join(B, 'clade.json')))[0]
p4.read(os.path.join(B, 'aa.nex'))
d = p4.Data(list(p4.var.alignments))
t = p4.var.trees[-1]
t.data = d
c0 = t.newComp(spec='specified', val=cfg['pi0'])
t.setModelComponentOnNode(c0, clade=1)
c1 = t.newComp(spec='specified', val=cfg['pit'])
h1 = [n for n in t.nodes if getattr(n, 'name', None) == 'H1'][0]
t.setModelComponentOnNode(c1, node=h1, clade=1)
t.newRMatrix(spec='lg')
t.setNGammaCat(nGammaCat=4)
t.newGdasrv(val=0.5)
t.setPInvar(val=0.0)
t.calcLogLike(verbose=0)
print('RES', repr(t.logLike))
