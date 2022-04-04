# tviblindi
This package is under development and depends on several libraries - issues during installation should be expected.

We recommend to pull the docker container provided with all dependencies and Rstudio server

# Docker container

port=7777 #e.g.\
data_path="path to data folder to mount"\
rpassword="password for rstudio server (user=rstudio)"\
docker run -it -d -p $port:8787 --name tviblindi_container -v $data_path:/data -e PASSWORD=$rpassword stuchly/tviblindi:v0422

then in your webrowser visit (user: rstudio, password: rpassword\
localhost:7777

# Installation 
devtools::install_github("stuchly/tviblindi")

Current version of tviblindi depends on CGAL library version 4.14

Macos:

brew tap-new CGAL/legacy   

brew extract --version=4.14 CGAL CGAL/legacy

brew install CGAL/legacy/CGAL@4.14  

# Citation
The article is in process

The shiny app was developed as a part of master's thesis of co-author [David Novak](https://github.com/davnovak):  [Studying lymphocyte development using mass cytometry](https://dspace.cuni.cz/handle/20.500.11956/119793?locale-attribute=en)
