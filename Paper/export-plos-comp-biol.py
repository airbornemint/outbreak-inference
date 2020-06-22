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
for idx, m in enumerate(re.finditer(r'\\paperincludegraphic\{\\figuresDir/(.*)\}', tex)):
	os.rename(f'{ m.group(1) }.tif', f'Fig{ idx + 1 }.tif')
