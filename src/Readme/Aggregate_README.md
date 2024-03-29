[Back to main README](../../README.md#tab-documentation)

[//]: # (document start)

## Aggregate
Since we only show a subset of experiments in the scatter plot, we provide an additional visualization that gives an aggregated view of the data in the back-end. 

![aggregate screenshot](https://user-images.githubusercontent.com/45741696/227915323-e7c7a913-3ace-4f1a-b539-2e3aaed27479.PNG)



### Select Feature (orange)
Users can choose a feature they want to visualize using the auto-complete text-input field. The values of the chosen feature are then shown in the background of the scatter plot as an aggregation of the data in the back-end into hexagonal bins. The aggregated value for each hexagon is encoded with colors.

Depending on whether or not the chosen feature is of temporal nature (i.e., temporal features are features that have values for each experiment cycle), a slider is shown, where users can choose the timestep they want to see. 

### Color encoding (yellow)
Users can adapt the color encoding by selecting a different colormap in the drop-down selection input.
Depending on whether the chosen feature consists of one or two values, the color encoding is slightly different:

![aggregate screenshot](https://user-images.githubusercontent.com/45741696/227915961-06f467d8-8ae6-4d00-b0be-62afae0e23b5.PNG)

![aggregate screenshot](https://user-images.githubusercontent.com/45741696/227916038-0a0db7bc-4f2d-47a9-aae8-0c46666ccb10.gif)


		
For features consisting of one value, the feature values are linearly mapped to the chosen color scale. 

If a feature consists of two variables (e.g., mean and variance of a prediction) we use a bi-variate mapping or the [VSUP (Value-Suppressing Uncertainty Palettes)](http://idl.cs.washington.edu/papers/uncertainty-palettes/) mapping. The former mapping linearly encodes one value as hue and the second value as saturation. The encoding can be viewed as a rectangular grid showing the hue on one axis and the saturation on the other axis. VSUP considers one of the values as a measure of “uncertainty”. With increasing uncertainty, the separation of values becomes less important and therefore doesn’t need as many distinguishable colors as those values that have a higher certainty (for more information, check out the paper that introduced [VSUP](http://idl.cs.washington.edu/papers/uncertainty-palettes/)). This mapping of the two values can be viewed as a “wedge” shape with the value encoded as hue and the uncertainty encoded as saturation.

### Settings (green)
Users can adjust the settings of the aggregation visualization. By default, the range used for mapping values to colors is derived from the data (i.e., the minimum and maximum value of the selected feature). Users can also customize the range, which can be useful, for example, if users want to compare several timesteps or features with each other.

![aggregate screenshot](https://user-images.githubusercontent.com/45741696/227916404-017014c4-2e2a-4e51-8647-682092211da8.PNG)

Finally, users can also choose the aggregation function that should be used to calculate the aggregated value. 

![aggregate screenshot](https://user-images.githubusercontent.com/45741696/227916514-9e2dac6c-6de2-4f01-8b96-1ba338113c0b.PNG)

If the selected feature has two values users can also choose, which of the two values should be considered as uncertainty for the color encoding. 

