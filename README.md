# ChemInformatics Model Explorer for Reaction Optimization (CIME4R): TODO…
This is the repository for CIME4R (as discussed in the paper). It builds upon the [Projection Space Explorer library](https://github.com/jku-vds-lab/projection-space-explorer).

#### Published in: TODO ####
#### DOI: TODO ###

This repository includes:
* The implementation of CIME4R
    * [Front-end](Application/) web application written in TypeScript using React
    * [Back-end](Application/reaction_cime/) python server
* [Documentation](#documentation)
* [Installation](#installation)
* [How to cite?](#how-to-cite)

Check out our [paper](TODO url) for further details about the implementation and use cases of CIME4R. 

Check out the [DEMO website](https://reaction-optimization.jku-vds-lab.at) of CIME4R, which includes the datasets used in the use cases.

Check out the [dataset generation examples](TODO add examples) if you want to try CIME4R with your own dataset.

Check out the [example datasets](TODO upload to osf) from the paper's use cases.

TODO: add table of contents

# Documentation CIME4R
The ChemInformatics Model Explorer for Reaction Optimization (short CIME4R) extension of the [Projection Space Explorer library](https://github.com/jku-vds-lab/projection-space-explorer/tree/develop) allows users to interactively explore the parameter space of chemical reactions and information about the iterative optimization process. The application allows users to understand how a machine learning model arrives at its decision on which experiments to perform next in retrospect (e.g., as proposed in [EDBO](https://www.nature.com/articles/s41586-021-03213-y)). It also facilitates interactive human-AI collaboration for reaction optimization to combine the advantages of both worlds for final decision-making: AI precision and human/expert intuition.
With CIME4R, users can apply a 2D projection to the provided reaction optimization data and show the high-dimensional data in a LineUp table, parallel coordinates plot, or hover view.
Since parameter spaces of chemical reactions can be huge, users can apply filtering or random subsampling to only show a (representative) subset of the data. The remaining data can be shown optionally using an aggregated view of the projected data. 
Users can interactively select data points (note: each data point represents one experiment configuration and will be called “experiment” in this documentation) in a 2D scatter plot and show summary statistics of features of all selected experiments in a summary visualization. 
Instructions for installing the application are provided at the end of this documentation.

## General/Controls
This section explains the general layout of the tool and the basic controls with which you can interact with the tool.

### View Components
- Left Menu Drawer (orange): Shows tabs that contain different groups of actions
- Projection (2-dimensional) View (yellow): Shows the current projection of the data and allows the user to interact with the low-dimensional projection of the experiments
- Tabular (high-dimensional) View (green): Can be dragged up from the bottom of the window to show a [LineUp](https://lineup.js.org/) table or parallel coordinates plot of the high dimensional space of the experiments

<img src="https://user-images.githubusercontent.com/45741696/227914093-57b73b24-308a-4dc8-afca-2668972fb42a.PNG" width="700">

### Controls
The following describes a list of controls:
- hover over item: shows a detailed view of the item
- left-click on item: select this item
- left-click + strg on item: toggle the selection status (i.e. if the item is selected, it is removed from selection; if the item is not selected, it is added to the selection)
- left-click on group-center: select the whole group
- left-click + strg on group-center: add the group to the selection
- left-click + drag on group-center: draw a directed edge to another group center
- left-click + drag: new lasso selection of experiments
- left-click + strg + drag: toggles the selection (i.e. unselected points that are within the lasso are added to the selection and selected points that are within the lasso are deselected)
- right-click + drag: allows you to move the whole scatter plot
- mouse wheel: zoom in and out to get a more/less detailed view of the experiments in the scatter plot
- right-click on background or item: opens a context menu with various options
- right-click on group center: opens group context menu that allows users to delete a group or start the group comparison feature
- right-click on group edge: opens group context menu that allows deleting the edge


## Dataset
In the “Dataset” tab users can choose to either upload their own dataset (orange) or load datasets that were already uploaded previously (yellow).
The “Select dataset” lists datasets that are already available in the back-end (from any user!) and can be deleted with the delete button next to the filename.
The list can also be manually refreshed with the refresh button next to “Select dataset” (this is only necessary if another user uploads a file during a simultaneous session and the current user needs this exact file).
If a user wants to upload a custom file, it must be a CSV or a zipped CSV. CIME4R recognizes special naming of columns in the dataset as described in the “Data Format” subsection. We provide a [datafile generation example](TODO: add examples) to get users started with their own datasets.
Finally, in the advanced settings (green) users can specify a SMILES lookup table and export/import a PSE session. 
With the SMILES lookup table, users can define key-value pairs of SMILES strings and human-readable names for those SMILES, which are then used by CIME4R to show the human-readable names instead of SMILES. The file has to be a CSV file with the columns “smiles” and “shortname”. Check out the [example SMILES lookup table](TODO: add example).
The export button allows users to save the current session of CIME4R so that they can later continue or even share the session with a colleague.


### Data Format
Data is handed to the system using a [comma separated values (CSV)](https://en.wikipedia.org/wiki/Comma-separated_values) file format that contains a collection of chemical reaction parameter configurations, measured target values (e.g., measured yield), during which cycle the measurement was performed, and additional properties that can be customized (e.g., predicted yield, SHAP values).
New files are first uploaded to the python back-end that runs with Flask (https://palletsprojects.com/p/flask/) and then preprocessed and stored in a [PostgreSQL](https://www.postgresql.org/) database.
For big files, the initial upload and preprocessing can take several minutes. If the files are already uploaded, it is much faster.

An example dataset can be found in [TODO: add example dataset](TODO). Datasets used in CIME4R's article are available in the data repository: [TODO: add link to OSF](TODO)

### Special Column Names
Users can define arbitrary column names. There are some special column names and column modifiers that are recognized by the system and have special meanings:
- Specifying the columns **x** and **y** tells the system to initialize the scatter plot according to these values. Otherwise, the x and y coordinates are randomly initialized and can later be overwritten with projections.
- The column **groupLabel** specifies the group each experiment belongs to. 
- The column **experiment_cycle** specifies the cycle in which this experiment setting was measured. If this experiment was not yet performed, the value should be **-1**.
 
### Column Modifiers
Column modifiers allow the system to give columns semantic meaning. The following gives a list of modifiers and how to use them:
- “measured”: this modifier indicates the target column. When the target is **yield**, the column should be specified as **measured_yield**. 
- “exp_param”: modifier that indicates that the column is an experiment parameter (e.g., "sulfonyl_equiv_exp_param", "temperature_exp_param")
- “smiles”: modifier that indicates that the column contains a SMILES string. SMILES are then automatically represented as a 2d structure (e.g., “sulfonyl_fluoride_SMILES_index”)
- “desc”: this modifier specifies that the column contains a descriptor value for a chemical compound. Usually, a descriptor consists of a lot of values (i.e., columns) that create overhead that is not important to explicitly show to users (i.e., they should not be shown in the table view), but are nice to use for projection.
- time-series data: since we deal with cyclic data (i.e., experiments are done in an iterative process), we can specify new values of a cycle as additional columns. The system automatically detects time-series data when they end with an **underscore** followed by a **number** (e.g., test_0, test_1, test_2). All columns that have the same name without the number are recognized to belong together. The number indicates the order.
- “pred”, “predicted”: are special modifiers for time-series data that indicate that these values were predicted by a model. These modifiers tell the system that this time-series data contains tuples. The naming must start with "pred" or “predicted” followed by the name of the variable and end with the timestep number, all separated with underscores (e.g. pred_mean_0, pred_var_0,...)
- “shap”: another special modifier for time-series data. It tells the system that this column contains diverging values as created by [SHAP (SHapley Additive exPlanations)](https://shap.readthedocs.io/en/latest/). The column names must end with "_shap" followed by the timestep (e.g., temperature_shap_0).


## Projection
When the data is loaded the x and y columns are used as initial positions for the scatter plot. If x and y are not specified they will be randomly initialized. 
The values for x and y can then be re-calculated with a projection method. 

Currently, there are options for [UMAP](https://umap-learn.readthedocs.io/), [t-SNE](https://opentsne.readthedocs.io/), and [PCA](https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html) projection implemented. The projection calculations are done in the back-end on the whole dataset. The positions of the points can be refined with [Overlap removal](https://arxiv.org/pdf/1903.06262.pdf) ([Code on Github](https://github.com/fpaulovich/dimensionality-reduction)), which uses a grid-based method to reduce overlaps of points and therefore helps to minimize visual clutter in the scatter plot.


### Parameters (orange)
For calculating the projection, users can click on one of the projection methods (orange). This opens a dialog, where users can choose the features which should be used for the projection by selecting and deselecting the corresponding checkboxes. To select or deselect whole semantic groups of features (e.g., select all features that belong to a descriptor, or all experiment parameters), users can interact with the checkboxes next to the group name. Users can also collapse and expand the list of experiments in a group.

The range value indicates the minimum and maximum values of the feature.

Users can also choose whether or not a numerical feature should be normalized. For most distance metrics, normalization applies standardization (i.e., subtract by mean and divide by standard deviation). Currently only when choosing [Gower’s distance](https://github.com/wwwjk366/gower), the values are normalized to the range [0; 1]. 

The “Weight” column can be used to change the weight for each feature, which determines the impact of a feature on the projection. For now, the weighting option is only considered in combination with Gower’s distance. By default, all features have equal weights (i.e., 1). Users can either change the weighting for each feature individually, or they can set the weights for an entire group of features. If the weight value for the entire group is changed, this value is evenly distributed among all features in this group.

The distance metric is automatically adapted when changing the selected features. If the selected features only include numerical columns, the default metric is “Euclidean”. If the columns are categorical only, “Jaccard” distance is chosen. If the columns contain mixed data types we use “Gower” distance.

Furthermore, users can adjust other hyperparameters used for the projection. Noteworthy here is the checkbox **Seed Position**, which tells the system to initialize the projection with the current positions of the experiments instead of using a random initialization.

Parameters that cannot be defined by the user are set to the defaults suggested in the corresponding frameworks (https://opentsne.readthedocs.io/en/latest/, https://umap-learn.readthedocs.io/en/latest/, https://scikit-learn.org/stable/modules/generated/sklearn.decomposition.PCA.html).  

### Progress (yellow)
The “Projection” tab panel includes a view that shows the progress of a projection as soon as the projection starts to calculate. Here, the calculations can be canceled by the user. 
Projections are calculated in the back-end for all available experiments. The points shown in the scatter plot automatically update their positions, as soon as the calculations are finished.

### Settings (green)
Apart from using the calculated projection coordinates, users can also choose features for the x and y positions of the points. By default, however, the projection coordinates are used (i.e., “x” and “y).


## Filter
As already mentioned previously the points shown in the scatter plot are only a subset of the whole dataset. This is necessary because reaction datasets can become very large, which is a challenge for computational resources on the one hand, and human processing capabilities on the other hand (i.e., showing all points to a user can become overwhelming quickly). 

The subset available in the front-end can be adjusted by users in the *Filter* tab panel. 

### Filter Info (orange)
The filter info tells users how many experiments are shown in the front-end and the total number of points in the back-end. It also indicates this information as a percentage value.

### Filter Settings (yellow)
By default, the filter is set to show those experiments for which we have actual measurement values. Whether or not measurement values are available for an experiment is indicated by the **experiment_cycle** feature. Therefore, the filter for this feature is set to only include experiments that have a positive value for the experiment_cycle. 

Users can adjust this filter by 
**adding** a new feature: using the auto-complete input field to select a feature
**adjusting** the value of a feature that should be filtered: for numerical features, users can set the minimum and maximum range of a feature on a slider; for categorical features, users can (de-)select discrete values of the feature
**removing** a feature from the filter: using the delete button next to the feature; this feature will not be considered for filtering anymore.

When the desired setting is found, users have to click on “Apply filter”, only then the data in the front-end will be updated and the filter settings will be stored in the back-end. Sometimes users might set the filters too loose resulting in too many experiments shown in the front-end. We defined an upper limit of 10,000 points that can be loaded into the front-end. If this limit is exceeded, a random subsample of 10,000 points is loaded to the front-end. Also, a hint is shown that tells users that only a subset of the filtered points can be seen. If users are not happy with their current settings, they can use the “Reset filter”. This will reset the filter settings to default settings in the back-end.

During the exploration of experiments, users might cross a specific experiment that is of interest to them. If they want to look into experiments with the same parameters as the found experiment, they can simply right-click on that point and choose the “Set filters to feature values of this experiment” context item. This automatically sets all filters to show exactly this experiment. Users can then adjust the feature values to their interests by removing certain constraints. This function aids users to arrive at the desired filter settings more quickly.

### Exception settings (green)
Sometimes users might be interested in a certain region of the scatter plot and would like to see all experiments. In that case, **Exceptions** can be defined. They are defined by the x and y coordinates and a range around which all experiments are loaded to the front-end. Exceptions are shown despite the filter settings (e.g., the filter is set to only show experiments that have measured values, except for the defined region in which all experiments are shown).

Users can define such exceptions with a right-click on the scatter plot. This opens the context menu showing the item “Show experiments around this position”. Clicking on this context item adds this position to the exception list. This function allows users to get more detailed information about an area that is of special interest to them.


## Aggregate
Since we only show a subset of experiments in the scatter plot, we provide an additional visualization that gives an aggregated view of the data in the back-end. 


### Select feature (orange)
Users can choose a feature they want to visualize using the auto-complete text-input field. The values of the chosen feature are then shown in the background of the scatter plot as an aggregation of the data in the back-end into hexagonal bins. The aggregated value for each hexagon is encoded with colors.

Depending on whether or not the chosen feature is of temporal nature (i.e., temporal features are features that have values for each experiment cycle), a slider is shown, where users can choose the timestep they want to see. 

### Color encoding (yellow)
Users can adapt the color encoding by selecting a different colormap in the drop-down selection input.
Depending on whether the chosen feature consists of one or two values, the color encoding is slightly different:
		
For features consisting of one value, the feature values are linearly mapped to the chosen color scale. 

If a feature consists of two variables (e.g., mean and variance of a prediction) we use a bi-variate mapping or the [VSUP (Value-Suppressing Uncertainty Palettes)](http://idl.cs.washington.edu/papers/uncertainty-palettes/) mapping. The former mapping linearly encodes one value as hue and the second value as saturation. The encoding can be viewed as a rectangular grid showing the hue on one axis and the saturation on the other axis. VSUP considers one of the values as a measure of “uncertainty”. With increasing uncertainty, the separation of values becomes less important and therefore doesn’t need as many distinguishable colors as those values that have a higher certainty (for more information, check out the paper that introduced [VSUP](http://idl.cs.washington.edu/papers/uncertainty-palettes/)). This mapping of the two values can be viewed as a “wedge” shape with the value encoded as hue and the uncertainty encoded as saturation.

### Settings (green)
Users can adjust the settings of the aggregation visualization. By default, the range used for mapping values to colors is derived from the data (i.e., the minimum and maximum value of the selected feature). Users can also customize the range, which can be useful, for example, if users want to compare several timesteps or features with each other.


Finally, users can also choose the aggregation function that should be used to calculate the aggregated value. 

If the selected feature has two values users can also choose, which of the two values should be considered as uncertainty for the color encoding. 



## Encoding
In the "Encoding" tab panel users can change the marks and channels of the displayed data.
- **shape by**: select a categorical feature with less than six distinct values and encode each value as a different mark 
- **opacity by**: select a numerical feature and scale the opacity of each point by that value; the upper and lower limit of the brightness can be adjusted with the scale below; if nothing is selected, the slider can be adjusted to set the general opacity value of all points
- **size by**: select a numerical feature and scale the size of each point by that value; the upper and lower limit of the size can be adjusted with the scale below; if nothing is selected, the slider can be adjusted to set the general size value of all points
- **color by**: select a categorical or numerical feature that defines the color of the points; the colormap can be chosen below and depends on whether the feature is numerical or categorical


## Selection Info
In this tab panel summary visualizations of selected points are shown. The info box at the top (orange) shows the number of selected experiments and the total experiments in the front-end.

Below this information, there are three icon buttons (yellow) explained from left to right:
External summary: open the selection visualization in an external window
Clear selection: reset the selection to no selected points
Summary configuration: opens a dialog with settings where users can (i) choose features to show in the visualization, and (ii) where to show the hover view (either in the top-right or bottom-left corner of the screen)

### Selection Visualization (green)
This visualization shows summary plots for a set of chosen features. 
For numerical features, the visualization is a density plot that gives an overview of the value distribution of selected experiments (blue) and all experiments in the scatter plot (black outline).
For categorical features, a bar chart gives an overview of the counts for each discrete value of the selected experiments (blue) and all experiments in the scatter plot (black outline). The bars are sorted by the count of the selected experiments.
Users can hover over the density plots or bar charts to see further details.
The plots are sorted by their purity. In other words, the homogeneity of feature values in the subset of selected experiments is calculated and used for sorting.

When selecting a hexagon (showing the aggregated data in the back-end), the plots use summary data from the back-end (i.e., blue areas correspond to the experiments within the selected hexagon and black outlines represent all experiments in the dataset).

### Hover View
The hover view is similar to the selection visualization; it shows a summary of the features of the experiment that is currently hovered in the scatter plot.


## Groups
In the "Groups" tab panel users can adjust group settings, automatically define groups by clustering and choose between group collections.



### Group Settings (orange)
A toggle allows users to show or hide experiments in the scatter plot. The other toggle allows users to show or hide group centers (grey diamonds).

If users click on a group center (grey diamond), all experiments belonging to that group are highlighted. Depending on which option users chose for the “Group visualization” the item membership of the selected group is visualized differently.  If **Contour plot** is selected, the experiments belonging to that group are surrounded by contour lines. If **Star visualization** is selected, there are lines drawn from the group center to each item. If **None** is selected, the points belonging to the group are just highlighted.

### Define groups by clustering (yellow)
Automatic Clustering of the projected features can be done in this panel. The algorithm used for clustering is [HDBSCAN](https://hdbscan.readthedocs.io/). 
Parameters can be changed either by adjusting the slider (few groups...many groups), or by enabling the **Advanced**-Mode. Chosen parameters are always synchronized with the values in the advanced user inputs. Any other possible parameters that could be used for HDBSCAN are set to the default parameters that can be retrieved from the HDBSCAN docs.


### Groups and collections (green)
A group collection is a set of groups - and possible connections between those groups - that were either created automatically or manually composed. This way, users can view various groupings by just switching between collections.

A new group collection can be created by clicking **New**. 
Users can manually add groups to a new or existing collection by selecting points in the scatter plot and choosing "Define group from selection" from the context menu that opens with a right-click on the scatter plot.
The group collection can be renamed and deleted in the “Settings” dialog.

The groups in a collection are listed below the user select. Each item in the list represents one group. If a user clicks on a group, the corresponding points are highlighted in the scatter plot.
Holding CTRL adds a group to the selection.
Next to each group label there is a settings button where users can adjust group names or delete a group.

### Group sequences and comparison
Users can manually define relationships between groups by drawing directed edges between group centers in the scatter plot (orange). Edges can be removed again with a right-click on the edge and selecting “Delete edge”. 


When the desired relationship tree/graph is built, users can use the “Group comparison” to clearly show differences between connected groups. 
To activate the “Group comparison” view, users must right-click on one of the group centers to open the context menu (yellow). From the context menu users can either choose “Group sequences … starting from this group” or “Group sequences … between 2 groups”. The first option shows all paths and comparisons starting from the chosen group. The latter option only shows the path between the two selected groups.

In the “Group comparison” view, the group summaries (for details see [Selection visualization](TODO: add link to this section)) are shown on the left part of the view and the comparison between two groups is shown on the right part of the view (green).

The differences for numerical features are encoded with two box plots. The upper box plot shows the feature value distribution of the first group and the lower box plot shows the distribution of the second group. 
Differences for categorical values are encoded with bar charts. Each categorical value that changes its occurrence count between the two groups is visualized as one bar that either goes in a positive or negative direction. If the count of this feature value is lower in the second group, the bar shows the percentage of the difference in the negative direction. If the count is higher in the second group, the bar shows the percentage of the difference in the positive direction.
The difference plots are sorted by how much a feature changes between the two groups.

## Tabular view
For high-dimensional data exploration, we included a [LineUp table](https://lineup.js.org/) and Parallel Coordinates Plot, which can be viewed on-demand. 

### LineUp
To show the table you can either use one of the buttons in the tab panel, or you can drag the component from the bottom of the window to increase the size of the table. 


The table shows all properties that were included in the provided dataset except properties that have the "desc" modifier (i.e., descriptor features are usually too numerous to show). 

For time series data, CIME4R combines all steps belonging to the same feature. LineUp shows the time series data in one column per feature in the form of a heatmap (the exact visualization type can be chosen as usual in LineUp).

For time series data consisting of two values (e.g., the mean and variance of a prediction) we implemented a special column that shows a line chart that draws one value as the line over time and the other value as an area around this line (i.e., the line shows the mean prediction, which is surrounded by a confidence interval using the variance).

Using the "smiles" modifier, users can manually specify, which properties represent SMILES strings. For each column that contains SMILES, the column shows the 2D structure of the SMILES feature.
The SMILES columns have some additional features:
- Users can filter those columns by substructure (a valid SMILES string must be provided in the filter input).
- Changing the width of those columns dynamically adapts row heights, which provides a better view of the 2D structures.
- When grouping several rows, this column displays the maximum common substructure of all compounds in the group.

All LineUp functionalities are included like filtering, searching, sorting, etc.
The grouping functionality can be performed in all columns, especially relevant is group by selected experiments and group by group labels, which actively uses features of the Projection Space Explorer.

The table can be used interactively with the scatter plot that represents the projected space and the summary view that shows selected experiments:
- Hovering experiments in the table highlights the corresponding experiments in the other views as well and vice versa.
- Users can select experiments in the table, which are also selected in the other views and vice versa.



### LineUp Settings
The **View all experiments** button automatically makes the table component visible - if it was not shown yet - and removes all filters. 

The **View selected experiments** button automatically makes the table component visible - if it was not shown yet - and filters the table by the selected experiments. 

The **Show cell values** toggle can be enabled to show values in numerical table cells. If it is disabled, the values are only shown for highlighted rows.

### Parallel Coordinates
Users can switch to a parallel coordinates plot by selecting **Parallel Coordinates** in the tab panel.

The usual functionalities of parallel coordinates plots are available, like defining and adjusting range constraints on the axis, or rearrangement of axis positions.

The parallel coordinates view is linked to the scatter plot: choosing constraints in parallel coordinates also highlights the corresponding points in the scatter plot and vice versa.

### Parallel Coordinates Settings
The **Choose attributes** button lets users choose, which features to show in the parallel coordinates plot.

The **Reset constraints** button removes all constraints set by the user.

Finally, users can export and import constraints for later use.


# Installation

There are multiple ways to run CIME4R. Option 1 is the easiest method.

## Option 1 (recommended) - Run CIME4R with Docker
Once you have Docker installed, you can quickly run the following commands and have CIME4R ready to use.

To **install** the latest version of CIME4R: 
TODO: ask michael about this
TODO: add docker for database
```bash 
docker pull ghcr.io/jku-vds-lab/reaction-cime:develop
docker run -d -p 9000:9000 --name cime4r --detach jkuvdslab/cime
```

To **update** CIME4R:
```bash
docker rm --force cime
docker pull jkuvdslab/cime
docker run -d -p 8080:8080 --name cime --detach jkuvdslab/cime
```

To **uninstall** CIME4R:
```bash
docker rm --force cime
```


## Option 2 - Build and run CIME4R from source (for development)

### Back-end
First, create a new virtual environment for the dependencies
```bash
python -m venv .venv
```
and activate it
```bash

# Ubuntu
source .venv/bin/activate

# Windows (cmd)
.\.venv\Scripts\activate
```

Then install all dependencies (including dev dependencies)

```bash
make develop
```

start the database with docker
```
docker compose up
```

and finally start the python server
```bash
python reaction_cime
```

### Front-end
Make sure you have the latest yarn version installed (`corepack enable` when using Node 16).

First install the required dependencies with yarn
```bash
yarn install
```

and launch the webpack-dev-server via
```bash
yarn start
```

Now, if a login screen pops up, you can use admin:admin to login. If you want to disable the login screen and go directly to the application, create a `reaction_cime/.env` with the following contents. After restarting the server, you will be automatically logged in.

```
VISYN_CORE__SECURITY__STORE__NO_SECURITY_STORE__ENABLE=true
VISYN_CORE__SECURITY__STORE__NO_SECURITY_STORE__USER=admin
```

### Link PSE
If you want to make changes to PSE and view the changes without having to push to the repo and reinstalling dependencies, the recommended way is to use the yarn link/portal and/or our webpack resolveAliases feature.

First, clone `projection-space-explorer` into the current directory (i.e. into the reaction-cime directory). Do not install `projection-space-explorer`, as we do not want any `node_modules` within that folder, as it should use the ones from the reaction-cime directory.

Add a portal to the local `projection-space-explorer` in the reaction-cime package.json (**this is a local change and should not be committed!**):

```json
  "resolutions": {
    ...
    "projection-space-explorer": "portal:./projection-space-explorer"
  },
```

Now, install everything via `yarn install`.

To now include `projection-space-explorer` to your webpack build, add a `.yo-rc-workspace.json` and update the `resolveAliases` (note the single `.`, i.e. using the current folder):

```json
{
  "resolveAliases": {
    "projection-space-explorer": "./projection-space-explorer/src/index.ts"
  }
}
```

With that, you can now edit all files of `projection-space-explorer`, including auto-completion (as the node_modules of the application will be used as main lookup), and get hot-reloading.

## Option 3 - Run Application with Docker from Source
```bash
yarn install
```

```bash
yarn run webpack:prod
```

```bash
docker build -f Dockerfile -t reaction_cime .
```

```bash
docker run --rm -it --network host reaction_cime
```

Beware that you will need a Postgres to run the image. By default, it will use the connection string in `settings.py`, which you can override via ENV variables. For example, you can set `REACTION_CIME__DBURL=postgresql://...` and use any database of your liking.



# How to cite?

You can cite CIME4R using the following bibtex:

```bibtex
@article{humer2023cime4r,
  author={TODO},
  journal={TODO},
  title={TODO},
  year={2023},
  doi={TODO},
  volume={TODO},
  number={TODO},
}
```



