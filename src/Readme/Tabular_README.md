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
