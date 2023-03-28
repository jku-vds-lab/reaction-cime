[Back to main README](../../README.md#tab-documentation)

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
