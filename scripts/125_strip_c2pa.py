#!/usr/bin/env python3
"""125: remove C2PA provenance metadata that a file-transfer tool adds to figures on the way to the server (a caBX chunk in PNG files;
an xmlns:c2pa attribute and a <metadata><c2pa:manifest> block in SVG files), so the committed figures are byte-identical to what was rendered.
Usage: python3 scripts/125_strip_c2pa.py file.png file.svg ...   (edits in place; anything else in the file is left unchanged)"""
import re, struct, sys
def png(p):
    b = open(p, 'rb').read(); assert b[:8] == b'\x89PNG\r\n\x1a\n', p; out = b[:8]; i = 8
    while i < len(b):
        n = struct.unpack('>I', b[i:i+4])[0]; t = b[i+4:i+8]; ch = b[i:i+12+n]; i += 12 + n
        if t != b'caBX': out += ch
    open(p, 'wb').write(out)
def svg(p):
    s = open(p, encoding='utf-8').read()
    s = s.replace(' xmlns:c2pa="http://c2pa.org/manifest"', '')
    s = re.sub(r'<metadata>\s*<c2pa:manifest>.*?</c2pa:manifest>\s*</metadata>', '', s, flags=re.S)
    open(p, 'w', encoding='utf-8').write(s)
for p in sys.argv[1:]:
    png(p) if p.lower().endswith('.png') else svg(p); print('stripped', p)
