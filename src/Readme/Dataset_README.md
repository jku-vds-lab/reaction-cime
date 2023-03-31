[Back to main README](../../README.md#tab-documentation)

[//]: # (document start)

## Dataset
In the “Dataset” tab users can choose to either upload their own dataset (orange) or load datasets that were already uploaded previously (yellow).
The “Select dataset” lists datasets that are already available in the back-end (from any user!) and can be deleted with the delete button next to the filename.
The list can also be manually refreshed with the refresh button next to “Select dataset” (this is only necessary if another user uploads a file during a simultaneous session and the current user needs this exact file).
If a user wants to upload a custom file, it must be a CSV or a zipped CSV. CIME4R recognizes special naming of columns in the dataset as described in the “Data Format” subsection. We provide a [datafile generation example](TODO: add examples) to get users started with their own datasets.
Finally, in the advanced settings (green) users can specify a SMILES lookup table and export/import a PSE session. 
With the SMILES lookup table, users can define key-value pairs of SMILES strings and human-readable names for those SMILES, which are then used by CIME4R to show the human-readable names instead of SMILES. The file has to be a CSV file with the columns “smiles” and “shortname”. Check out the [example SMILES lookup table](TODO: add example).
The export button allows users to save the current session of CIME4R so that they can later continue or even share the session with a colleague.

![dataset screenshot](https://user-images.githubusercontent.com/45741696/227914723-6b7a48a0-9d41-4519-b4b6-7cbb31f4f525.PNG)

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