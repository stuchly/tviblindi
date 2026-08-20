# Linux test container (RStudio Server + Shiny GUI)

A native-Linux image for testing `tviblindi`'s Shiny GUI outside macOS - built
while tracking down a GUI bug that only reproduced on Linux, not macOS. This
Dockerfile is also what `stuchly/tviblindi:arm64` (see main `tviblindi/README.md`)
is published from - the main README's Docker section is what most users
want; this file is for building it yourself (e.g. to iterate on the image)
or for local dev/testing rather than pulling the published tag.

The original Docker image documented in the main README (`stuchly/tviblindi`,
no tag) explicitly does **not** support Apple Silicon ("Mac compatibility is
limited to Intel"); this one does, by building natively for `linux/arm64`
rather than relying on a prebuilt amd64 image under emulation.

## Build & run

From the `tviblindi/` package root (the build context needs the whole package
source):

```bash
docker build --platform linux/arm64 -f docker/linux-test/Dockerfile -t tviblindi-linux-test .
docker run -d --platform linux/arm64 --name tviblindi-linux-test \
  -p 8787:8787 -p 3838:3838 -e PASSWORD=<your-password> tviblindi-linux-test
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
`--platform linux/amd64` explicitly.

## Why this image needed its own set of fixes

None of these are needed for the conda-based build (`recipe/`,
`env/build_env.yaml`) - conda-forge's pinned `cgal`/`r-bh` packages already
avoid them. They only bit because this image installs R packages from CRAN
(via Posit Package Manager) and system libraries via `apt`, which float to
whatever's current on Ubuntu:

1. **`BH` must be pinned to `1.87.0-1`.** Newer `BH` bundles a Boost that
   dropped `boost::mpl::if_c`, which CGAL 5.x's headers (`number_utils.h`,
   `Sqrt_extension_type.h`, etc.) still reference - exact same issue as
   `PACKAGING_NOTES.md` #6 documents for the conda recipe, just hit via a
   different package manager. This is a CGAL-5.x-only concern: verified CGAL
   6.x compiles fine even against `BH` 1.90+, since it dropped that Boost
   pattern internally. `tviblindi` itself has no CGAL upper-version pin -
   apt's stock `libcgal-dev` (5.6 on Ubuntu Noble) is used as-is, no
   vendoring needed.
2. **The `CGAL_HEADER_ONLY` preprocessor macro must be defined when
   compiling tviblindi's own C++.** Apt's `libcgal-dev` ships no compiled
   `libCGAL.so` (verified via `dpkg -L libcgal-dev` - header-only by design),
   so `CGAL::Random` and a few other normally-non-inline CGAL symbols need
   this macro to become inline instead of requiring that missing library.
   Setting it via a `CPPFLAGS` environment variable does **not** work - R's
   build system on this image doesn't propagate that env var into the actual
   compile line (verified empirically). It's injected directly into
   `/usr/include/CGAL/config.h` instead.
3. **`libglpk-dev`** is needed for `igraph` to load at all (`libglpk.so.40`
   missing otherwise) - not a `tviblindi`-specific issue, just an apt package
   that isn't pulled in by `libcgal-dev`/friends.
4. **`libwebpmux3`** is needed for the R package `ragg` to load at all
   (pulled in transitively via `tidyverse`, which `tviblindi` Depends on).
   Without it, any Shiny plot fails with `unable to load shared object
   '.../ragg/libs/ragg.so': libwebpmux.so.3: cannot open shared object
   file` - found by actually launching the Shiny GUI and clicking into the
   "Terminal nodes selection" / "Homology classes by persistence selection"
   panels, not by any earlier compile or install step (`R CMD INSTALL`
   succeeds fine without it - `ragg` only gets touched at plot-render time).

## `rocker/rstudio` vs. a conda-based image + RStudio Server

An earlier attempt installed RStudio Server standalone on top of a
conda-based image (matching `env/build_env.yaml`). That's a dead end for
Apple Silicon: RStudio Server's official Linux `.deb` from Posit is
**amd64-only**, forcing the whole image under Rosetta emulation - which
turned out to be genuinely slow (`library(tviblindi)` took ~16s per fresh R
process, vs. 0.9s native). `rocker/rstudio` ships its own working arm64
RStudio Server build and handles the `PASSWORD` env var / session startup
correctly out of the box, so it was used here instead, at the cost of the
CGAL/BH fixes above (since it uses a system R + apt/CRAN packages, not a
conda-pinned CGAL/Boost pair).
