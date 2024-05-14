
This repository is associated with the article Zhang et al. (2021) "Meeting the Bare Minimum: Quality Assessment of Idiographic Temporal Networks Using Power Analysis and Predictive Accuracy Analysis".

To reproduce the simulation presented in the article, you need to follow the steps below:

1. **run "src/estimate_parameters.R":** it estimates/gets and exports the parameters values of the VAR(1) in the "data/df_param.csv" and "data/list_param.rda" files. Note that the data of Epskamp, van Borkulo, et al. (2018) are not provided in this repository.
2. **run "src/main.R":** run the power analysis and the PAA. export the results of the simulation in the "data/data_sim/" folder. All the custom functions used in the "main.R" script are in the "src/FUN" folder. 
3. **run data/preprocess_results.R:** merge the csv files stored in the "data/data_sim/" folder.
4. **knit "figures/create_fig.Rmd":** all the plots and tables of the results from the simulation are also exported in the "figures/outputs" folder.


# Authors

* Jordan Revol
* Yong Zhang

# Citation

The preprint associated with this repository is available at: xxx