<img src="vignettes/tviblindi_logo.png" alt="tviblindi" width=220>

**Topological and Geometrical Tools for Single-Cell Data**

*tviblindi* is a trajectory inference toolkit for single-cell data.

This package is under development and depends on several libraries - issues during installation are expected.

We recommend to pull the Docker container provided with all dependencies and an RStudio server.

- - - - -

Parts of code from following libraries are used:\
[RcppAnnoy](https://cran.r-project.org/web/packages/RcppAnnoy/index.html)\
[Gudhi](https://gudhi.inria.fr)\
[Phat](https://www.sciencedirect.com/science/article/pii/S0747717116300098)\
[CGAL](https://www.cgal.org)

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

To be able to use default dimensionality reduction, install *vaevictis*

```
reticulate::py_install("git+https://github.com/stuchly/vaevictis.git",pip=TRUE)
```

# Usage

We include sample code below to run the *tviblindi* pipeline on synthetic data.

```
library(tviblindi)
data(tviblindi_dyntoydata)
group_id<-tviblindi_dyntoydata[,1]
datainput<-as.matrix(tviblindi_dyntoydata[,-1])
tv1<-tviblindi(data=datainput,labels=group_id)
DimRed(tv1)
DimRed(tv1,method="umap")

Set_origin(tv1,label = "M4",origin_name = "M4_hitting_time")
Set_origin(tv1,label = "M4",origin_name = "M4_hitting_distance")
KNN(tv1)
Cluster(tv1,K=225) #kmeans clustering
Filtration(tv1) #default setting is very conservative, less simplices could be
# created with same resolution (e.g. Filtration(tv1,alpha2=1))

Pseudotime(tv1,weighted = FALSE,origin_name = "M4_hitting_time")
Walks(tv1,N=1000,origin_name = "M4_hitting_time")

Pseudotime(tv1,weighted = TRUE,origin_name = "M4_hitting_distance")
Walks(tv1,N=1000,origin_name = "M4_hitting_distance")
launch_shiny(tv1)
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
[Preprint](https://www.biorxiv.org/content/10.1101/2023.07.13.547329v2) describing and validating the method.

See `supplemetary_note.pdf` in vignette for technical background.

The Shiny app (GUI) was developed as part of [David Novak](https://github.com/davnovak)'s master thesis: [Studying lymphocyte development using mass cytometry](https://dspace.cuni.cz/handle/20.500.11956/119793?locale-attribute=en).
