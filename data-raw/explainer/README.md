# Explainer figures

Figures of the article `vignettes/articles/how-crywolf-works.Rmd`, written to
`vignettes/articles/figures/how-*.png`. Run from the package root:

```sh
Rscript data-raw/explainer/make_explainer.R
```

- `make_explainer.R` draws the figures of data with 'ggplot2' and 'patchwork', in the palette
  and theme of the plots of the package.
- `model.dot` and `results.dot` are the diagrams, drawn with Graphviz (`dot`). Function names are
  in IBM Plex Mono.

Values in the figures are illustrative. The P(occupied) values are computed in the script with
ψ = 0.5, p11 = 0.6 and p10 = 0.05. The mixture uses p11 = 0.25, so that the two groups overlap.
