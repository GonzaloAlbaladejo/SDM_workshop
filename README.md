# Species distribution models (SDMs) in R
## Applications and limitations
This code provides a quick and easy hands-on introduction to species distribution and niche modelling. The code contains all the different steps to perform these analyses:
  1. Data gathering and preparation
  3. Species data acquisition, processing and clean-up
  4. Exploratory niche analysis using the environmental niche factor analysis (ENFA)
  5. A diversity of supervised and unsupervised statistical methods for the analysis of species distributions and habitat suitability
  6. SDM model averaging and forecasting

Most of the code is presented both in its raw annotated format ".R" and as an executable R-markdown format "Rmd". Independently of the format the code works the same. A `Data` folder is provided along with the code and some extra information needed for the analysis. Similarly, a list of custom functions is included in the `Functions` folder. These functions allow for the automation of several aspects of SDM modelling and are designed to facilitate the execution of the code.
The code and data are all self-contained, which means that once you download the repository you need to keep the same structure of folders for the code to find the right routes to the data.
The execution order is marked by the numbers at the start of each R or Rmd file (e.g. first 0. Data_preparation.R, followed by 1. Process_Spp_Data.R, etc...)
The functionality of some functions is linked to the latest versions of their dependent packages, so if something is not working try to upload the packages. If the problem persists try to annotate and post it in this GitHub repository. Most of the base functions are frequently checked and managed for error control.
If you use any of this code or functions in some of your projects, please don't forget to reference this repository.
