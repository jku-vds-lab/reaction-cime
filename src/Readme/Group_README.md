[Back to main README](../../README.md#tab-documentation)

[//]: # (document start)

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
