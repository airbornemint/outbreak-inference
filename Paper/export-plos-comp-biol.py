# Helper for exporting for PLOS Comp Biol 
# - Inlines bibliography
# - Renames images to match the text order

import re
import sys
import os

# Inline the bibliography
tex = open(sys.argv[1]).read()
tex = re.sub(r'\\bibliography\{(.*)\}', lambda f: open(f'{f.group(1)}.bbl').read(), tex)
open(sys.argv[1], 'w').write(tex)

# Get the figure names (in order) and rename the images
figRE = r'(?s)\\begin\{figure\}.*?\\end\{figure\}'
subFigRE = r'(?s)\\begin\{subfigure\}.*?\\end\{subfigure\}'
graphicRE = r'\\paperincludegraphic\{\\figuresDir/(.*?)\}'

for figIdx, figM in enumerate(re.finditer(figRE, tex)):
	if re.search(subFigRE, figM.group(0)):
		for subFigIdx, subFigM in enumerate(re.finditer(subFigRE, figM.group(0))):
			m = re.search(graphicRE, subFigM.group(0))
			os.rename(f'{ m.group(1) }.tif', f'Fig{ figIdx + 1 }{ chr(subFigIdx + ord("a")) }.tif')
	else:
		m = re.search(graphicRE, figM.group(0))
		os.rename(f'{ m.group(1) }.tif', f'Fig{ figIdx + 1 }.tif')
