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

## Why this image needs a few extra fixes

None of these are needed for the conda-based build (`recipe/`,
`env/build_env.yaml`) - conda-forge's pinned `cgal`/`r-bh` packages already
avoid them. They only apply here because this image installs R packages
from CRAN (via Posit Package Manager) and system libraries via `apt`, which
float to whatever's current on Ubuntu:

1. **`BH` is pinned to `1.87.0-1`.** Newer `BH` bundles a Boost that dropped
   `boost::mpl::if_c`, which CGAL 5.x's headers (`number_utils.h`,
   `Sqrt_extension_type.h`, etc.) still reference - same issue
   `PACKAGING_NOTES.md` #6 documents for the conda recipe, hit here via a
   different package manager. This only matters for CGAL 5.x: CGAL 6.x
   compiles fine even against `BH` 1.90+. `tviblindi` itself has no CGAL
   upper-version pin - apt's stock `libcgal-dev` (5.6 on Ubuntu Noble) is
   used as-is, no vendoring needed.
2. **The `CGAL_HEADER_ONLY` preprocessor macro is defined** when compiling
   `tviblindi`'s own C++. Apt's `libcgal-dev` ships no compiled
   `libCGAL.so` (header-only by design), so `CGAL::Random` and a few other
   normally-non-inline CGAL symbols need this macro to become inline
   instead of requiring that missing library. A `CPPFLAGS` environment
   variable doesn't work for this - R's build system here doesn't propagate
   it into the actual compile line - so it's injected directly into
   `/usr/include/CGAL/config.h` instead.
3. **`libglpk-dev`** is needed for `igraph` to load at all (`libglpk.so.40`
   missing otherwise) - just an apt package that isn't pulled in by
   `libcgal-dev`/friends.
4. **`libwebpmux3`** is needed for the R package `ragg` (pulled in
   transitively via `tidyverse`) to load. Without it, Shiny plots fail with
   `unable to load shared object '.../ragg/libs/ragg.so': libwebpmux.so.3:
   cannot open shared object file` - this only surfaces at plot-render
   time, not during `R CMD INSTALL`.

## Why `rocker/rstudio` instead of RStudio Server on a conda-based image

RStudio Server's official Linux `.deb` from Posit is amd64-only, so
installing it on top of a conda-based image (matching
`env/build_env.yaml`) only works under Rosetta emulation on Apple Silicon -
which is genuinely slow (`library(tviblindi)` took ~16s per fresh R
process under emulation, vs. 0.9s native). `rocker/rstudio` ships its own
working arm64 RStudio Server build and handles the `PASSWORD` env var /
session startup correctly out of the box, at the cost of the CGAL/BH fixes
above (since it uses a system R + apt/CRAN packages, not a conda-pinned
CGAL/Boost pair).
