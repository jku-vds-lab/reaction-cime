[Back to main README](../../README.md#tab-documentation)

[//]: # (document start)

## Filter
As already mentioned previously the points shown in the scatter plot are only a subset of the whole dataset. This is necessary because reaction datasets can become very large, which is a challenge for computational resources on the one hand, and human processing capabilities on the other hand (i.e., showing all points to a user can become overwhelming quickly). 

The subset available in the front-end can be adjusted by users in the *Filter* tab panel. 

![filter screenshot](https://user-images.githubusercontent.com/45741696/227915177-39b14697-6d42-41dd-b8c1-d936d7637e86.PNG)


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
