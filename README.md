# Species distribution models (SDMs) in R
## Code description
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

## Scripts actions
### 0. Data_preparation/ PrepareR_env
This code downloads the climatic environmental data from {WorldClim], the geographic information, and any other additional information needed for the description and delimitation of the study area.

### 1. Process_Spp_data
Downloads the spatial records for the selected species from the digital repositories of [GBIF](https://www.gbif.org/) as well as the taxonomic information from [Itis](https://itis.gov/) and the [IUCN red-list of threatened species](https://www.iucnredlist.org/). This code also retrieves the range data for the selected species. This information has been previously downloaded from the IUCN Red-List [spatial data repositories]() and placed on the `Data/Spatial/Processed/Vector` folder. If you want to execute this code for other species you will need to download this information.

### 2. SDM_models
This scripts combines the environmental and species spatial records, the study area and the species ranges to run a series of SDM algoritmns:
	a. Envairomental Niche Factor Analysis (ENFA)
	b. Binomial Generalized Linear Models (GLM)
	c. 
In addition to this, this code prepares the set of pseudo-absence and background data following different spatial sampling strategies. Once processed the models are exported into an R-data (rds) object for further analysis and processing into the `/Results/SDM_models` folder.

### 3. Model averaging and forecasting
This script takes the previously processed models along with future climatic data to run species distribution predictions under different scenarios of climate change.



If you use any of this code or functions in some of your projects, please don't forget to reference this repository.
