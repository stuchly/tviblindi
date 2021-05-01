# tviblindi

**Single-cell trajectory inference toolkit**

- - - - -

`tviblindi` puts concepts from graph theory and algebraic topology to use for trajectory inference (TI) in high-dimensional biological data (cytometry, scRNA-seq, CITE-seq).

We provide easy-to-use tools for identifying potential developmental trajectories and grouping them in a classification tree.
This includes a graphical user interface that enables the user to

* inspect trajectories in different classes by topological relationships,

* view the trajectories drawn on a 2-dimensional layout of input data,

* track progression of expression levels for markers of interest at different stages of development,

* check the composition of tracked populations in terms of labelled cell populations and

* export enhanced FCS files for viewing results of the TI analysis in FlowJo or other gating software.

## Installation

`tviblindi` is distributed as an open-source `R` package, you can install it via `devtools`.
It includes compiled code integrated via `Rcpp`, and requires `OpenMP` for parallelisation.
The `CGAL` library version 4 is required (`tviblindi` has not been adapted to work with newer versions of `CGAL` yet), which is why we include `CGAL` in the `src` directory.

## Usage

You can run the `tviblindi` trajectory-inference pipeline easily using only two functions from `tviblindi`: the first for generating a `tviblindi` object and simulating tentative developmental pathways, and the second for opening an interactive graphical user interface that helps you discover and classify pathways. For simplicity, you can use a sample mass-cytometry human thymus dataset included with the package. All parameters of the analysis are described in the documentation for each function.

```
library(tviblindi)
data('Thymus_expression') # load thymus expression data (coordinate matrix)
data('Thymus_labels')     # load labels per cell
tv <- InferTrajectories(
  expression = Thymus_expression,
  labels_per_event = Thymus_labels,
  origin_label = 'CD34+CD1a- DN T precursor'
)
Interactive(tv)
```

Alternatively, you can run the constituent functions that comprise the trajectory-inference pipeline one-by-one.

```
tv <- tviblindi(data = Thymus_expression, labels = Thymus_labels)
SetOrigin(tv, 'CD34+CD1a- DN T precursor')
ConstructkNNG(tv)
Denoise(tv)
AddLayout(tv, method = 'umap')
Plot(tv)
Cluster(tv)
Filter(tv)
ComputePseudotime(tv)
SimulateRandomWalks(tv)
Interactive(tv)
```

To inspect the connectedness of different populations in your dataset based on Louvain clustering, use the 'connectome' functionality.

```
Connectome(tv)
```

Additionally, to make a part of your point cloud denser, you can sample extra data points from an artificial distribution simulated over the neighbourhood of existing points.
You can choose, for instance, to inflate a single specific population of interest, changing the result of simulating developmental pathways as a result.

```
Upsample(tv, idcs = 'CD4+CD3+CD1a+ SP')
```

On the other hand, you can use `Downsample` to reduce the number of points.
However, it is preferrable to use the full data whenever possible.

In the process of running your TI pipeline (after generating a 2-d layout), you can call `Plot` to view your data.

Each user-level function is documented, so finding out about all its parameters (or references to literature which describe relevant concepts), use `?FunctionName` to view the documentation.

## How `tviblindi` works

For an explanation of the inner workings of `tviblindi`, check out `supplementary.pdf`.

