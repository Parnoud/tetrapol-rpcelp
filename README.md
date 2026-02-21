# rp-celp (Fork)

**A fork of [TenGbps/rp-celp](https://github.com/TenGbps/rp-celp)**

This repository is a fork of the original [rp-celp](https://github.com/TenGbps/rp-celp) project. The purpose of this fork is to try make it work and use it on modern system.

## About the Original Project

[rp-celp](https://github.com/TenGbps/rp-celp) is a audio codec for Tetrapol communication . The original repository is maintained by [TenGbps](https://github.com/TenGbps).


## Installation

Orignal codec from brmlab team

```bash
git clone https://github.com/yourusername/rp-celp.git
cd rp-celp
python3 decoder.py my-files.json audio.wav
```
--> No audible

2nd way 

```bash
g++ -o rpcelp rpcelp.cpp
cat /tmp/tch_channel.json | grep '"type": "VOICE"' | grep -oE '"value": "[0-9a-f]{30}"' | cut -d \" -f 4 | ./blhexbit.py | sed -re "s/(.)(.)(.)(.) (.)(.)(.)(.) /\8\7\6\5\4\3\2\1/g" | sed -re "s/^([01]{20})/\1/g" | ./rpcelp > out.raw
ffmpeg -f s16le -ar 8000 -ac 1 -i out.raw ffmpeg.wav
```

--> audible

3 way (live decoding from mono channel)
```bash
g++ -o rpcelp rpcelp.cpp
tetrapol_dump -r 42000 -t TCH 2>/dev/null  g| grep '"type": "VOICE"' | grep -oE '"value": "[0-9a-f]{30}"' | cut -d \" -f 4 | ./blhexbit.py | sed -re "s/(.)(.)(.)(.) (.)(.)(.)(.) /\8\7\6\5\4\3\2\1/g" | sed -re "s/^([01]{20})/\1/g" | ./rpcelp | ffplay -f s16le -ar 8000 -i pipe:0
```

--> live decoding (note optimize, cpu goes HOT)