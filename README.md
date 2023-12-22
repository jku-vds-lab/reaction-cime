# CIME4R: Exploring iterative, AI-guided chemical reaction optimization campaigns in their parameter space
This is the repository for CIME4R. It builds upon the [Projection Space Explorer library](https://github.com/jku-vds-lab/projection-space-explorer).

#### Submitted to: Journal of Cheminformatics ####
#### Preprint: https://doi.org/10.26434/chemrxiv-2023-218lq ####
#### DOI: 10.26434/chemrxiv-2023-218lq ###

This repository includes:
* The implementation of CIME4R
    * Front-end web application written in TypeScript using React
    * [Back-end](reaction_cime/) python server
* [Documentation](#documentation-cime4r)
* [Installation](#installation)
* [How to cite?](#how-to-cite)

Check out our [paper](https://doi.org/10.26434/chemrxiv-2023-218lq) for further details about the implementation and use cases of CIME4R. 

Check out the [DEMO website](https://reaction-optimization.jku-vds-lab.at) of CIME4R, which includes the datasets used in the use cases.

Check out the [example datasets and generation scripts](https://osf.io/vda72/) used in the paper.

# Documentation CIME4R
The ChemInformatics Model Explorer for Reaction Optimization (short CIME4R) extension of the [Projection Space Explorer library](https://github.com/jku-vds-lab/projection-space-explorer/tree/develop) allows users to interactively explore the parameter space of chemical reactions and information about the iterative optimization process. The application allows users to understand how a machine learning model arrives at its decision on which experiments to perform next in retrospect (e.g., as proposed in [EDBO](https://www.nature.com/articles/s41586-021-03213-y)). It also facilitates interactive human-AI collaboration for reaction optimization to combine the advantages of both worlds for final decision-making: AI precision and human/expert intuition.
With CIME4R, users can apply a 2D projection to the provided reaction optimization data and show the high-dimensional data in a LineUp table, parallel coordinates plot, or hover view.
Since parameter spaces of chemical reactions can be huge, users can apply filtering or random subsampling to only show a (representative) subset of the data. The remaining data can be shown optionally using an aggregated view of the projected data. 
Users can interactively select data points (note: each data point represents one experiment configuration and will be called “experiment” in this documentation) in a 2D scatter plot and show summary statistics of features of all selected experiments in a summary visualization. 
Instructions for installing the application are provided at the end of this documentation.

##### Table of Contents  
[General/Controls](#generalcontrols)  
[Dataset](#dataset)  
[Projection](#projection)  
[Filter](#filter)  
[Aggregate](#aggregate)  
[Encoding](#encoding)  
[Selection Info](#selection-info)  
[Groups](#groups)  
[Tabular view](#tabular-view)

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


## Tab Documentation
[Dataset](src/Readme/Dataset_README.md)

[Projection documentation](src/Readme/Projection_README.md)

[Encoding documentation](src/Readme/Encoding_README.md)

[Group documentation](src/Readme/Group_README.md)

[Selection documentation](src/Readme/Selection_README.md)

[Tabular documentation](src/Readme/Tabular_README.md)

[Filter documentation](src/Readme/Filter_README.md)

[Aggregate documentation](src/Readme/Aggregate_README.md)

# Installation

There are multiple ways to run CIME4R. Option 1 is the easiest method.

## Option 1 (recommended) - Run CIME4R with Docker Compose
Once you have Docker and Docker Compose installed, you can quickly run the following commands and have CIME4R ready to use.

Simply copy the `docker-compose-demo.yml` from this repository (or checkout this repository) and run: 

```bash
docker compose -f ./docker-compose-demo.yml up
```

This will start the server and a corresponding postgres for the datasets. You can the navigate to http://localhost:9000/ and use the application.

## Option 2 - Run CIME4R with Docker
If you don't want to use Docker Compose, you need to bring your own postgres. Make sure it is installed locally or override it with the env variable (REACTION_CIME__DBURL=postgresql://admin:admin@db_postgres:5432/db).

To **install** the latest version of CIME4R: 
```bash 
docker pull ghcr.io/jku-vds-lab/reaction-cime:develop
docker run -d -p 9000:9000 --name cime4r --detach jkuvdslab/cime
```


## Option 3 - Build and run CIME4R from source (for development)

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
  author={Humer, Christina and Nicholls, Rachel and Heberle, Henry and Heckmann, Moritz and Pühringer, Michael and Wolf, Thomas and Lübbesmeyer, Maximilian and Heinrich, Julian and Hillenbrand, Julius and Volpin, Giulio and Streit, Marc},
  journal={ChemRxiv},
  title={{CIME4R}: {Exploring} iterative, {AI}-guided chemical reaction optimization campaigns in their parameter space},
  shorttitle = {{CIME4R}},
  year={2023},
  month = dec,
  doi={10.26434/chemrxiv-2023-218lq},
  note={This content is a preprint and has not been peer-reviewed.}
}
```



