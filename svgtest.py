#!/usr/bin/env python3
# -*- coding: iso-8859-15 -*-

import svgwrite

d = svgwrite.Drawing(filename='images/alignment.svg',
                     size=(60,60))
circle = d.circle((30,30), 30, fill='red')
text = d.text('Test', (30,30),
              style='text-anchor:middle;\
                     dominant-baseline:middle',
              font_size='17px')

d.add(circle)
d.add(text)
d.save()
