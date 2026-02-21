#!/usr/bin/python3

import sys, binascii

for l in sys.stdin:
  for c in l:
    if c in "0123456789ABCDEFabcdef":
      d = binascii.unhexlify("0"+c)
      bits = bin(d[0])[2:]
      while len(bits) < 4:
        bits = "0"+bits
      print(bits, end=" ")
  print()
