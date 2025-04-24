#!/usr/bin/env python
import argparse

parser = argparse.ArgumentParser(
    description="Generate initial conditions for Taylor–Green vortex")
parser.add_argument("-l", "--level", type=int, required=True)
parser.add_argument("-o", "--output", type=str, required=True)
args = parser.parse_args()

import math
import numpy as np

L = 2 * math.pi
nvars = 4
n = 1 << args.level
dx = L / n
a = np.memmap(args.output, dtype=float, mode="w+", shape=(nvars, n, n, n))
U, V, W, P = a
i = np.arange(n)
s = np.sin(dx * i)
c = np.cos(dx * i)
np.einsum('i,j,k',  s, c, c, out=U)
np.einsum('i,j,k', -c, s, c, out=V)
P.fill(0)
W.fill(0)
