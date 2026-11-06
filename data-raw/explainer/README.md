# Explainer figures

Panels of the article `vignettes/articles/how-crywolf-works.Rmd`, written to
`vignettes/articles/figures/how-*.png`:

```sh
data-raw/explainer/make_explainer.sh
```

- `panels.py` draws each panel as an SVG string (`svg_helpers.py` holds the drawing helpers).
- `build.py` puts them in `explainer.html`, one figure per panel, with `base.css` and the
  Figtree and IBM Plex Mono fonts (downloaded from npm on the first run).
- `render.js` saves each figure as a PNG at 2x with Playwright.

Values in the drawings are illustrative; the P(occupied) values of the first two panels are
computed with ψ = 0.5, p11 = 0.6 and p10 = 0.05.
