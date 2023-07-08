<img src="vignettes/tviblindi_logo.png" alt="tviblindi" width=220>

**Topological and Geometrical Tools for Single-Cell Data**

*tviblindi* is a trajectory inference toolkit for single-cell data.

This package is under development and depends on several libraries - issues during installation are expected.

We recommend to pull the Docker container provided with all dependencies and an RStudio server.

- - - - -

Parts of code from following libraries are used:\
[RcppAnnoy](https://cran.r-project.org/web/packages/RcppAnnoy/index.html)\
[Gudhi](https://gudhi.inria.fr)\
[Phat](https://www.sciencedirect.com/science/article/pii/S0747717116300098)

- - - - -

<kbd>
  <img src="vignettes/tviblindi_workflow.png">
</kbd>


`tviblindi` puts concepts from graph theory and algebraic topology to use for trajectory inference (TI) in high-dimensional biological data (cytometry, scRNA-seq, CITE-seq).

We provide easy-to-use tools for identifying potential developmental trajectories and grouping them in a classification tree.
This includes a graphical user interface that enables the user to

* inspect trajectories in different classes by topological relationships,

* view the trajectories drawn on a 2-dimensional layout of input data,

* track progression of expression levels for markers of interest at different stages of development,

* check the composition of tracked populations in terms of labelled cell populations and

* export enhanced FCS files for viewing results of the TI analysis in FlowJo or other gating software.

# Docker container

(Apple Silicon is not supported for the time being: Mac compatibility is limited to Intel.)

With Docker installed, run the following code in a Unix terminal.

```
port=7777\
data_path="path to data folder to mount"\
rpassword="password for rstudio server (user=rstudio)"\
docker run -it -d -p $port:8787 --name tviblindi_container -v $data_path:/data -e PASSWORD=$rpassword stuchly/tviblindi
```

Then navigate to `localhost:7777` in your web browser.
Enter the default credentials when prompted (user: `rstudio`, password: `rpassword`).

`localhost` may also depend on your Docker daemon setting.

# Direct installation

Currently, *tviblindi* depends on the *CGAL* library version 4.14 (not higher).

See installation instructions [here](https://doc.cgal.org/4.14/Manual/installation.html) or follow the instructions below on Intel Macs.

```
brew tap-new CGAL/legacy   
brew extract --version=4.14 CGAL CGAL/legacy
brew install CGAL/legacy/CGAL@4.14
```

To install *tviblindi* in R, run

```
devtools::install_github('stuchly/tviblindi')
```

# Usage

We include sample code below to run the *tviblindi* pipeline on synthetic data.
R package *TDA* is required to make the synthetic dataset.

```
library(tviblindi)
sn           <- make_snowman3d(2000) # produces synthetic 3-dimensional 'snowman'
colnames(sn) <- c('A', 'B', 'C')     # simulated markers
lab          <- rep(1, nrow(sn))
lab[which.max(sn[,2])] <- 0          # choose snowman's head as 'cell of origin'
lab          <- as.factor(lab)

tv1                <- tviblindi(data = sn, labels = lab)
tv1$origin$default <- which.max(sn[,2])

KNN(tv1, 50, method = 'balltree')    # create k-NNG (balltree faster for small data)
if (FALSE) {
    Denoise(tv1)                     # reduce noise before witness complex construction (for real world data)
}
Cluster(tv1, K=225)       # k-means clustering (15*15 clusters)
Filtration(tv1)                      # witness complex filtration
DimRed(tv1)                          # create lower-dimensional embedding
Pseudotime(tv1, sym = 'min')
Walks(tv1, N = 1000)

launch_shiny(tv1)                    # note the 'Min trajectory count % per leaf' slider (to see all classes) & 'Point size' selector
```

<kbd>
  <img src="vignettes/tviblindi_gui.png">
</kbd>

To inspect the connectedness of different populations in your dataset based on Louvain clustering, use the fully automated 'connectome' functionality: produces a PNG file.

```
Connectome(tv1)
```

<center>
<kbd>
  <img src="vignettes/connectome.png" width=350>
</kbd>
</center>

# Reference
This package is under construction and we're working on a paper describing the method and validating it.

See `supplemetary_note.pdf` in vignette for some technical background (will be updated).

The Shiny app (GUI) was developed as part of [David Novak](https://github.com/davnovak)'s master thesis: [Studying lymphocyte development using mass cytometry](https://dspace.cuni.cz/handle/20.500.11956/119793?locale-attribute=en).
