[Back to main README](../../README.md#tab-documentation)

[//]: # (document start)

## Encoding
In the "Encoding" tab panel users can change the marks and channels of the displayed data.
- **shape by**: select a categorical feature with less than six distinct values and encode each value as a different mark 
- **opacity by**: select a numerical feature and scale the opacity of each point by that value; the upper and lower limit of the brightness can be adjusted with the scale below; if nothing is selected, the slider can be adjusted to set the general opacity value of all points
- **size by**: select a numerical feature and scale the size of each point by that value; the upper and lower limit of the size can be adjusted with the scale below; if nothing is selected, the slider can be adjusted to set the general size value of all points
- **color by**: select a categorical or numerical feature that defines the color of the points; the colormap can be chosen below and depends on whether the feature is numerical or categorical


## Selection Info
In this tab panel summary visualizations of selected points are shown. The info box at the top (orange) shows the number of selected experiments and the total experiments in the front-end.

![encoding screenshot](https://user-images.githubusercontent.com/45741696/227916658-23fa58c1-579d-431a-9e4d-45ab7cade39e.PNG)


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
