# tviblindi
This package is under development and issues during installation should be expected.

We recommend to pull the docker container provided with all dependencies and Rstudio server

# Installation 
devtools::install_github("stuchly/tviblindi")

Current version of tviblindi depends on CGAL library version 4.14

Macos:

brew tap-new CGAL/legacy   

brew extract --version=4.14 CGAL CGAL/legacy

brew install CGAL/legacy/CGAL@4.14  

# Citation
The article is in process

The shiny app was develloped as part of Master's thesis of David Novak [https://github.com/davnovak]: https://dspace.cuni.cz/handle/20.500.11956/119793?locale-attribute=en
