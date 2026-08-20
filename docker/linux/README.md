# Linux Docker image (RStudio Server + Shiny GUI)

Builds a Linux image with `tviblindi`, RStudio Server, and `vaevictis` -
this is what `stuchly/tviblindi:arm64` (see the main `tviblindi/README.md`
"Docker container" section) is published from. Most users should just pull
that published tag; build from this Dockerfile yourself if you want to
customize the image or verify a change before it's published.

The other image documented in the main README (`stuchly/tviblindi`, no tag)
is Intel/amd64 only. This one builds natively for `linux/arm64` as well, so
it also covers Apple Silicon hosts without relying on amd64 emulation.

## Build & run

From the `tviblindi/` package root (the build context needs the whole package
source):

```bash
docker build --platform linux/arm64 -f docker/linux/Dockerfile -t tviblindi-linux .
docker run -d --platform linux/arm64 --name tviblindi-linux \
  -p 8787:8787 -p 3838:3838 -e PASSWORD=<your-password> tviblindi-linux
```

- RStudio Server: `http://localhost:8787` (user `rstudio`, password from `PASSWORD`)
- To also get the Shiny GUI running automatically with the bundled
  `tviblindi_dyntoydata` example (KNN → Cluster → Filtration → Pseudotime →
  Walks → UMAP layout), run `Rscript /start_app.R` inside the container (or
  from an RStudio terminal) - it's copied in but not run by default, so a
  plain RStudio session doesn't have to wait for it. Once running, the app is
  on `http://localhost:3838`.
- If you only need RStudio (no auto-launched Shiny app), you don't need to
  touch `/start_app.R` at all - just connect and `library(tviblindi)`
  yourself.

For a build targeting an Intel/amd64 host or CI runner, drop
`--platform linux/arm64` (defaults to the host's own architecture) or pass
`--platform linux/amd64` explicitly - nothing in the Dockerfile itself is
architecture-specific.

## What's in the image

- Base: `rocker/rstudio:4.5` (R 4.5 + RStudio Server)
- `tviblindi` and its R dependencies, installed from CRAN via Posit Package
  Manager binaries
- CGAL: Ubuntu's apt package (5.6 on Noble) - `tviblindi` has no upper CGAL
  version pin, this is just whatever Ubuntu currently ships
- `vaevictis`, pinned to a specific commit, in its own Python venv
  (`/opt/venv-tviblindi`), with TensorFlow

See the Dockerfile's own comments for why each of the less obvious lines
(the `BH` pin, `CGAL_HEADER_ONLY`, etc.) is there.
