[Back to main README](../../README.md#tab-documentation)

[//]: # (document start)

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
