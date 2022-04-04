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


