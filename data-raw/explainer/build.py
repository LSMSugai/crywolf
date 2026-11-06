"""Write explainer.html: one figure per panel of the article “How crywolf works”."""
import os, base64
from panels import problem, model, mixture, psite, designs, threshold, covs, flow

d = os.path.dirname(os.path.abspath(__file__)) + '/'
fp = d + 'fontpkg/node_modules/@fontsource/'
ff = ''
for fam, pkg, ws in [('Figtree', 'figtree', [400, 700]), ('IBM Plex Mono', 'ibm-plex-mono', [400])]:
    for w in ws:
        b = base64.b64encode(open(f'{fp}{pkg}/files/{pkg}-latin-{w}-normal.woff2', 'rb').read()).decode()
        ff += f"@font-face{{font-family:'{fam}';font-weight:{w};font-style:normal;font-display:block;src:url(data:font/woff2;base64,{b}) format('woff2')}}\n"
css = open(d + 'base.css').read()
own = '''
body{padding:16px;background:#fff}
.stage{--c:var(--s4);--soft:var(--s4-soft)}
figure.panel{margin:0 0 24px;display:block}
figure.w{width:820px} figure.h{width:500px}
'''
# id, function, title, wide
PANELS = [("problem", problem, "The problem", True),
          ("model", model, "The model: two stories for every site", False),
          ("mixture", mixture, "How the model tells them apart", False),
          ("p-occupied", psite, "How P(occupied) is computed for one site", True),
          ("designs", designs, "It needs some known truth: two designs", True),
          ("threshold", threshold, "Where the threshold enters", False),
          ("covariates", covs, "What can change p10", False),
          ("results", flow, "From detections to results", True)]
figs = "".join(f'<figure class="panel {"w" if w else "h"}" id="{pid}"><div class="phead"><span class="ttl">{t}</span></div><div>{f()}</div></figure>'
               for pid, f, t, w in PANELS)
html = (f'<!doctype html><html lang="en"><head><meta charset="utf-8"><title>crywolf explainer</title>'
        f'<style>{ff}{css}{own}</style></head><body><section class="stage st4">{figs}</section></body></html>')
open(d + 'explainer.html', 'w').write(html)
print('ok explainer:', ", ".join(p[0] for p in PANELS))
