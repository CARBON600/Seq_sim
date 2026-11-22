import json
import os
import sys

import p4

from _env import directory

W = directory('SEQ_SIM_EXT_WORK', 'a writable working directory')
B = os.path.join(W, 'b2')
cfg = json.load(open(os.path.join(B, 'hky.json')))[int(sys.argv[1])]
p4.read(os.path.join(B, 'dna.nex'))
d = p4.Data(list(p4.var.alignments))
t = p4.var.trees[-1]
t.data = d
t.newComp(spec='specified', val=cfg['pi'])
t.newRMatrix(spec='2p', val=cfg['kappa'])
t.setNGammaCat(nGammaCat=1)
t.setPInvar(val=0.0)
t.calcLogLike(verbose=0)
print('RES', cfg['tag'], repr(t.logLike))
