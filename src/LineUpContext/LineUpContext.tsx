import * as React from 'react';
import { connect, ConnectedProps } from 'react-redux';
import * as LineUpJS from 'lineupjs';
import './LineUpContext.scss';
import {
  IStringFilter,
  createSelectionDesc,
  Column,
  ERenderMode,
  IDynamicHeight,
  IGroupItem,
  Ranking,
  IRenderContext,
  IOrderedGroup,
  ICellRenderer,
  ICellRendererFactory,
  IDataRow,
  IGroupCellRenderer,
  renderMissingDOM,
  StringColumn,
  deriveColumnDescriptions,
} from 'lineupjs';
import * as _ from 'lodash';
import {
  AStorytelling,
  DiscreteMapping,
  EXCLUDED_COLUMNS,
  PrebuiltFeatures,
  selectVectors,
  setHoverState,
  AShallowSet,
  mapValueToColor,
} from 'projection-space-explorer';
import * as d3v5 from 'd3v5';
import isEqual from 'lodash.isequal';
import { arrayEquals, mapSmilesToShortname } from '../Utility/Utils';

import { ReactionCIMEBackendFromEnv } from '../Backend/ReactionCIMEBackend';
import { TestColumn } from './LineUpClasses/TestColumn';
import { setLineUpInputLineup } from '../State/LineUpInputDuck';
import { AppState } from '../State/Store';

/**
 * Declares a function which maps application state to component properties (by name)
 *
 * @param state The whole state of the application (contains a field for each duck!)
 */
const mapStateToProps = (state: AppState) => ({
  lineUpInput: state.lineUpInput,
  lineUpInput_data: state.dataset?.vectors,
  lineUpInput_columns: state.dataset?.columns,
  currentAggregation: state.currentAggregation,
  activeStory: AStorytelling.getActive(state.stories),
  pointColorScale: state.multiples.multiples.entities[state.multiples.multiples.ids[0]]?.attributes.pointColorScale,
  channelColor: state.multiples.multiples.entities[state.multiples.multiples.ids[0]]?.attributes.channelColor,
  detailView: state.detailView,
  // splitRef: state.splitRef
  // hoverState: state.hoverState
});

/**
 * Declares a function which maps dispatch events to component properties (by name)
 *
 * @param dispatch The generic dispatch function declared in redux
 */
const mapDispatchToProps = (dispatch) => ({
  setCurrentAggregation: (samples: number[]) => dispatch(selectVectors(samples)),
  setLineUpInput_lineup: (input) => dispatch(setLineUpInputLineup(input)),
  setHoverstate: (state, updater) => dispatch(setHoverState(state, updater)),
});

/**
 * Factory method which is declared here so we can get a static type in 'ConnectedProps'
 */
const connector = connect(mapStateToProps, mapDispatchToProps);

/**
 * Type that holds the props we declared above in mapStateToProps and mapDispatchToProps
 */
type PropsFromRedux = ConnectedProps<typeof connector>;

/**
 * Type that holds every property that is relevant to our component, that is the props declared above + our OWN component props
 */
// type Props = PropsFromRedux & {
//     onFilter: any
//     // My own property 1
//     // My own property 2
// }

type Props = PropsFromRedux & {
  // onFilter;
};

// let lineup = null;
const UPDATER = 'lineup';
export const UNIQUE_ID = 'unique_ID';

/**
 * Our component definition, by declaring our props with 'Props' we have static types for each of our property
 */
export const LineUpContext = connector(function ({
  lineUpInput,
  lineUpInput_data,
  lineUpInput_columns,
  currentAggregation,
  channelColor,
  setCurrentAggregation,
  setLineUpInput_lineup,
  // onFilter,
  activeStory,
  pointColorScale,
  setHoverstate,
  detailView,
}: Props) {
  // In case we have no input, dont render at all
  if (!lineUpInput || !lineUpInput_data || !detailView.open) {
    // splitRef?.current?.setSizes([100, 0])
    return null;
  }
  const lineupRef = React.useRef<any>();

  // eslint-disable-next-line react-hooks/exhaustive-deps
  const debouncedHighlight = React.useCallback(
    _.debounce((hover_item) => setHoverstate(hover_item, UPDATER), 200),
    [],
  );

  const preprocessLineupData = (data) => {
    // if (activeStory)
    //   ACluster.deriveVectorLabelsFromClusters(
    //     data,
    //     Object.values(activeStory.clusters.entities)
    //   );
    const lineupData = new Array<any>();
    const columns = {};
    data.forEach((element) => {
      // if(element[PrebuiltFeatures.ClusterLabel].length <= 0){
      //     element[PrebuiltFeatures.ClusterLabel] = [-1];
      // }
      const row = {};

      for (const i in lineUpInput_columns) {
        const col = lineUpInput_columns[i];

        if (!EXCLUDED_COLUMNS.includes(i) && (Object.keys(col.metaInformation).length <= 0 || !col.metaInformation.noLineUp)) {
          if (Object.keys(col.metaInformation).length > 0 && col.metaInformation.timeSeriesGroup) {
            // if(col.metaInformation.timeSeriesGroup.endsWith(""))
            const split = col.metaInformation.timeSeriesGroup.split(':');
            if (split.length <= 1) {
              // if the string is separated with a colon, only the first part of the string is considered as the group. the second part of the string determines a sub value of this group
              if (Object.keys(row).includes(col.metaInformation.timeSeriesGroup)) {
                row[col.metaInformation.timeSeriesGroup].push(element[i]);
              } else {
                row[col.metaInformation.timeSeriesGroup] = [element[i]];
                columns[col.metaInformation.timeSeriesGroup] = col;
                columns[col.metaInformation.timeSeriesGroup].metaInformation.listData = true;
                columns[col.metaInformation.timeSeriesGroup].metaInformation.range = col.metaInformation.globalRange; // TODO: iterate over all columns and derive global min/max
                columns[col.metaInformation.timeSeriesGroup].metaInformation.colorMapping = col.metaInformation.colorMapping;
              }
            } else {
              const group_name = split[0];
              const var_name = split[1];
              if (Object.keys(row).includes(group_name)) {
                if (Object.keys(row[group_name]).includes(var_name)) {
                  row[group_name][var_name].push(element[i]);
                } else {
                  row[group_name][var_name] = [element[i]];
                }
              } else {
                row[group_name] = {};
                row[group_name][var_name] = [element[i]];
              }

              // update column metaInformation
              if (Object.keys(columns).includes(group_name)) {
                columns[group_name].metaInformation.globalMin = Math.min(columns[group_name].metaInformation.globalMin, element[i]);
                columns[group_name].metaInformation.globalMax = Math.max(columns[group_name].metaInformation.globalMax, element[i]);
              } else {
                columns[group_name] = col;
                columns[group_name].metaInformation.customLineChart = true;
                columns[group_name].metaInformation.globalMin = element[i];
                columns[group_name].metaInformation.globalMax = element[i];
              }
            }
          } else if (Object.keys(col.metaInformation).length > 0 && col.metaInformation.lineup_meta_column) {
            // add meta column, that can be used by other columns
            row[i] = element[i];
            columns[i] = col;

            row[col.metaInformation.lineup_meta_column] = element[i];
            columns[col.metaInformation.lineup_meta_column] = { metaInformation: { hide: true } };
          } else {
            row[i] = element[i];
            columns[i] = col;
          }
        }
      }

      row[PrebuiltFeatures.ClusterLabel] = element[PrebuiltFeatures.ClusterLabel].toString();
      row[UNIQUE_ID] = element.__meta__.meshIndex;
      lineupData.push(row);

      // console.log(element)
      // let row = Object.assign({}, element)
      // row[PrebuiltFeatures.ClusterLabel] = element[PrebuiltFeatures.ClusterLabel].toString();
      // row[UNIQUE_ID] = element["__meta__"]["meshIndex"];
      // lineup_data.push(row);
    });
    return [lineupData, columns];
  };

  const clearAutomaticFilters = (lineUpInput, filter) => {
    if (filter) {
      for (const key in filter) {
        const { lineup } = lineUpInput;
        const ranking = lineup.data.getFirstRanking();
        if (key === 'selection') {
          const filter_col = ranking.children.find((x) => {
            return x.desc.column === UNIQUE_ID;
          });
          filter_col?.clearFilter();
        } else {
          const filter_col = ranking.children.find((x) => {
            return x.desc.column === key;
          });
          filter_col?.clearFilter();
        }
      }
    }
  };
  const getLineupDump = (lineUpInput) => {
    if (lineUpInput.lineup) {
      clearAutomaticFilters(lineUpInput, lineUpInput.filter);
      const dump = lineUpInput.lineup.dump();
      return dump;
    }
    return null;
  };

  React.useEffect(() => {
    // if(lineUpInput.dump){
    //     try {
    //         const json_parsed = JSON.parse(lineUpInput.dump)
    //         const restored = fromDumpFile(json_parsed)
    //         console.log(restored);
    //         const builder = buildLineup(lineUpInput.columns, restored.dat).restore(restored.dump);
    //         // const builder = LineUpJS.builder(restored.data).restore(restored.dump);
    //         lineup?.destroy();
    //         lineup = builder.build(lineup_ref.current);
    //         return;
    //     } catch (error) {
    //         console.log(error);
    //     }
    // }

    const tempData = preprocessLineupData(lineUpInput_data);
    const lineupData = tempData[0];
    const columns = tempData[1];

    const builder = buildLineup(columns, lineupData, pointColorScale, channelColor); // lineUpInput_data
    const dump = getLineupDump(lineUpInput);

    lineUpInput.lineup?.destroy();
    let lineup;
    lineup = builder.buildTaggle(lineupRef.current);
    if (dump) {
      lineup.restore(dump);
    }

    const ranking = lineup.data.getFirstRanking();

    // add selection checkbox column
    let selectionCol = ranking.children.find((x) => x.label === 'Selection Checkboxes');
    if (!selectionCol) {
      selectionCol = lineup.data.create(createSelectionDesc());
      if (selectionCol) {
        ranking.insert(selectionCol, 1);
      }
    }

    // // make lineup filter interact with the scatter plot view
    // ranking.on('orderChanged.custom', (previous, current, previousGroups, currentGroups, dirtyReason) => {

    //     if (dirtyReason.indexOf('filter') === -1) {
    //         return;
    //     }

    //     const onRankingChanged = (current) => {
    //         for (let i=0; i < lineUpInput.data.length; i++) {
    //             lineUpInput.data[i].view.lineUpFiltered = !current.includes(i);
    //         }

    //         onFilter()

    //     }

    //     onRankingChanged(current)
    // });

    // make lineup selection interact with the scatter plot view
    lineup.on('selectionChanged', (currentSelection_lineup) => {
      // if(currentSelection_lineup.length == 0) return; // selectionChanged is called during creation of lineup, before the current aggregation was set; therefore, it would always set the current aggregation to nothing because in the lineup table nothing was selected yet

      const currentSelectionCcatter = lineUpInput_data
        .map((x, i) => {
          if (x.__meta__.selected) return i;
          return undefined;
        })
        .filter((x) => x !== undefined);

      if (!arrayEquals(currentSelection_lineup, currentSelectionCcatter)) {
        // need to check, if the current lineup selection is already the current aggregation
        const agg = new Array<number>();
        currentSelection_lineup.forEach((index) => {
          agg.push(lineUpInput_data[index].__meta__.meshIndex);
        });

        setCurrentAggregation(agg);
      }
    });

    lineup.on('highlightChanged', (idx) => {
      let hoverItem;
      if (idx >= 0) {
        hoverItem = lineUpInput_data[idx];
      }
      debouncedHighlight(hoverItem);
    });

    // update lineup when smiles_column width changes
    if (smilesStructureColumns.length > 0 || customChartColumns.length > 0) {
      const custom_chart_cols = ranking.children.filter((x: any) => customChartColumns.includes(x.label));
      for (const i in custom_chart_cols) {
        const custom_chart_col = custom_chart_cols[i];
        custom_chart_col.on('widthChanged', (prev, current) => {
          lineup.update();
        });
      }

      const lineupSmilesCols = ranking.children.filter((x: any) => smilesStructureColumns.includes(x.label));
      for (const i in lineupSmilesCols) {
        const lineup_smiles_col = lineupSmilesCols[i];
        lineup_smiles_col.on('widthChanged', (prev, current) => {
          lineup.update();
        });

        // custom filter adapted from michael
        const filterChanged = (prev, cur: IStringFilter) => {
          if (prev?.filter !== cur?.filter) {
            // only update, if it is a new filter
            const filter = typeof cur?.filter === 'string' ? cur?.filter : null; // only allow string filters -> no regex (TODO: remove regex checkbox)
            if (lineup_smiles_col && filter) {
              ReactionCIMEBackendFromEnv.getSubstructureCount(
                lineUpInput_data.map((d) => d[lineup_smiles_col.desc.column]),
                filter,
              )
                .then((matches) => {
                  const validSmiles = matches.filter(([smiles, count]) => count > 0).map(([smiles, count]) => smiles);
                  lineup_smiles_col.setFilter({
                    filter,
                    valid: new Set(validSmiles),
                    filterMissing: cur.filterMissing,
                  });
                })
                .catch((e) => {
                  lineup_smiles_col.setFilter(null);
                });
            }
          }
        };
        lineup_smiles_col.on(StringColumn.EVENT_FILTER_CHANGED, filterChanged);
      }
    }

    setLineUpInput_lineup(lineup);

    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [lineUpInput_data, lineUpInput_columns, activeStory, activeStory?.clusters, activeStory?.clusters?.ids.length, lineUpInput.update]);

  // React.useEffect(() => { //TODO: not working...
  //     // update lineup, if current storybook (current cluster) changed
  //     if(lineUpInput.lineup){
  //         const data_provider = lineUpInput.lineup.data;
  //         let lineup_data = preprocess_lineup_data(lineUpInput_data);
  //         console.log("setdata")
  //         data_provider.setData(lineup_data);

  //         const ranking = lineUpInput.lineup.data.getFirstRanking();
  //         const my_col_builder = LineUpJS.buildCategoricalColumn(PrebuiltFeatures.ClusterLabel);
  //         console.log(lineup_data)
  //          ranking.insert(lineUpInput.lineup.data.create(my_col_builder.build(lineup_data)));
  //         ranking.insert(lineUpInput.lineup.data.create(my_col_builder.build([]])));

  //         // const ranking = lineUpInput.lineup.data.getFirstRanking();
  //         // // let cluster_col = ranking.columns.find(x => x.desc.column == PrebuiltFeatures.ClusterLabel);
  //         // // const my_desc = cluster_col.desc;
  //         // // my_desc.categories = ["test"]
  //         // // const my_col = new CategoricalColumn(cluster_col.id, cluster_col.desc)
  //         // const my_col_builder = LineUpJS.buildCategoricalColumn(PrebuiltFeatures.ClusterLabel);
  //         // // console.log()
  //         // ranking.insert(lineUpInput.lineup.data.create(my_col_builder.build(lineup_data))); //asSet(',')

  //         // // data_provider.setData(lineup_data)
  //         // // lineUpInput.lineup.update();
  //         // // lineUpInput.lineup?.setDataProvider(data_provider);
  //         // lineUpInput.lineup.restore(lineUpInput.lineup.dump())

  //         // // console.log(cluster_col.dump())
  //         // // console.log(lineUpInput.lineup.dump())

  //     }
  // }, [activeStory, activeStory?.clusters, activeStory?.clusters?.length]);

  // this effect is allways executed after the component is rendered when currentAggregation changed
  React.useEffect(() => {
    if (lineUpInput.lineup != null) {
      // select those instances that are also selected in the scatter plot view
      if (currentAggregation.aggregation && currentAggregation.aggregation.length > 0) {
        const currentSelectionCcatter = lineUpInput_data
          .map((x, i) => {
            if (x.__meta__.selected) return i;
            return undefined;
          })
          .filter((x) => x !== undefined);
        lineUpInput.lineup.setSelection(currentSelectionCcatter);

        // const lineup_idx = lineup.renderer?.rankings[0]?.findNearest(currentSelection_scatter);
        // lineup.renderer?.rankings[0]?.scrollIntoView(lineup_idx);

        // set the grouping to selection checkboxes -> uncomment if this should be automatically if something changes
        // const ranking = lineup.data.getFirstRanking();
        // let selection_col = ranking.children.find(x => x.label == "Selection Checkboxes");
        // ranking.groupBy(selection_col, -1) // remove grouping first
        // ranking.groupBy(selection_col);
      } else {
        lineUpInput.lineup.setSelection([]);
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [lineUpInput.lineup, currentAggregation]);

  React.useEffect(() => {
    if (lineUpInput.lineup && lineUpInput.lineup.data) {
      const ranking = lineUpInput.lineup.data.getFirstRanking();
      clearAutomaticFilters(lineUpInput, lineUpInput.previousfilter);
      if (lineUpInput.filter) {
        for (const key in lineUpInput.filter) {
          const curFilter = lineUpInput.filter[key];

          if (key === 'reset' && curFilter) {
            ranking.clearFilters();
          } else if (key === 'selection') {
            const filter_col = ranking.children.find((x) => {
              return x.desc.column === UNIQUE_ID;
            });

            let regexStr = '';
            lineUpInput.filter[key].forEach((element) => {
              regexStr += '|';
              regexStr += element; // ["__meta__"]["meshIndex"];
            });
            regexStr = regexStr.substr(1); // remove the leading "|"
            const myRegex = new RegExp(`^(${regexStr})$`, 'i'); // i modifier says that it's not case sensitive; ^ means start of string; $ means end of string
            filter_col?.setFilter({
              filter: myRegex,
              filterMissing: true,
            });
          } else {
            const filterCol = ranking.children.find((x) => {
              return x.desc.column === key;
            });
            const myRegex = new RegExp(`^(.+,)?${curFilter}(,.+)?$`, 'i'); // i modifier says that it's not case sensitive; ^ means start of string; $ means end of string
            filterCol?.setFilter({
              filter: myRegex,
              filterMissing: true,
            });
          }
        }
      }
    }
    // eslint-disable-next-line react-hooks/exhaustive-deps
  }, [lineUpInput.lineup, lineUpInput.filter]);

  // https://github.com/lineupjs/lineup_app/blob/master/src/export.ts
  return (
    <div className="LineUpParent">
      <div
        style={{
          clear: 'both',
          position: 'absolute',
          top: '1px',
          bottom: 0,
          left: 0,
          right: 0,
          padding: 0,
        }}
        ref={lineupRef}
        id="lineup_view"
      />
    </div>
  );
});

const WIDTH_HEIGHT_RATIO = 2;
let smilesStructureColumns = new Array<any>();
let customChartColumns = new Array<any>();
function myDynamicHeight(data: IGroupItem[], ranking: Ranking): IDynamicHeight {
  if (smilesStructureColumns.length > 0) {
    const cols = ranking.children.filter((x) => smilesStructureColumns.includes(x.label) || customChartColumns.includes(x.label));

    if (!cols || cols.length === 0) return { defaultHeight: 25, height: () => 25, padding: () => 0 };

    const colHeights = cols.map((x) => {
      if (customChartColumns.includes(x.label)) return x.getWidth() / WIDTH_HEIGHT_RATIO; // for chart, the width should be bigger than the height
      return x.getWidth(); // for images it is square
    });
    const colHeight = Math.max(Math.max(...colHeights), 25); // col.getWidth();

    const height = function (item: IGroupItem | Readonly<IOrderedGroup>): number {
      return colHeight;
    };
    const padding = function (item: IGroupItem | Readonly<IOrderedGroup>): number {
      return 0;
    };
    return { defaultHeight: colHeight, height, padding };
  }
  return { defaultHeight: 25, height: () => 25, padding: () => 0 };
}

// const base_color = undefined;
const baseColor = '#c1c1c1';
// const base_color = "#1f77b4";
function buildLineup(cols, data, pointColorScale, channelColor) {
  // console.log(channelColor) //TODO: update lineup colorscale, if sth changes; TODO: do this for all columns, not just groupLabel
  let groupLabelCatColor;
  if (channelColor?.key === PrebuiltFeatures.ClusterLabel) {
    // TODO: update colormapping code; does not work since colormapping of PSE changed...
    // let groupLabel_mapping = new DiscreteMapping(
    //   pointColorScale,
    //   new ShallowSet(
    //     data.map((vector) => vector[PrebuiltFeatures.ClusterLabel])
    //   )
    // );
    const groupLabelMapping = {
      scale: pointColorScale,
      values: AShallowSet.create(data.map((vector) => vector[PrebuiltFeatures.ClusterLabel])),
      type: 'categorical',
    } as DiscreteMapping;
    groupLabelCatColor = groupLabelMapping.values
      .filter((cat) => cat && cat !== '')
      .map((cat) => {
        return { name: cat, color: mapValueToColor(groupLabelMapping, cat).hex };
      });
  }

  const builder = LineUpJS.builder(data);

  for (const i in cols) {
    const col = cols[i];
    const show = true; //! (typeof col.metaInformation.hideLineUp !== 'undefined' && col.metaInformation.hideLineUp); // hide column if "hideLineUp" is specified -> there is a lineup bug with that option

    // if (!EXCLUDED_COLUMNS.includes(i) && (Object.keys(col.metaInformation).length <= 0 || !col.metaInformation.noLineUp)) { // only if there is a "noLineUp" modifier at this column or thix column is excluded, we don't do anything
    if (col.metaInformation.imgSmiles) {
      const smilesCol = `Structure: ${i}`;
      smilesStructureColumns.push(smilesCol);
      builder.column(
        LineUpJS.buildColumn('mySmilesStructureColumn', i)
          .label(smilesCol)
          .renderer('mySmilesStructureRenderer', 'mySmilesStructureRenderer')
          .width(50)
          .build([]),
      );
      // uncomment if you also want to show the smiles string, not just the structure
      // builder.column(LineUpJS.buildStringColumn(i).width(50).custom("visible", show).color(base_color));
    } else if (col.metaInformation.customLineChart) {
      // builder.column(LineUpJS.buildNumberColumn(i).label(i).asMap().renderer("myLineChartRenderer", "myLineChartRenderer").width(50).build([]));
      builder.column(
        LineUpJS.buildColumn('myLineChartColumn', i)
          .label(i)
          .custom('min', col.metaInformation.globalMin)
          .custom('max', col.metaInformation.globalMax)
          .renderer('myLineChartRenderer', 'myLineChartRenderer')
          .width(150)
          .build([]),
      );
      customChartColumns.push(i);
    } else if (i === PrebuiltFeatures.ClusterLabel) {
      const clustCol = LineUpJS.buildCategoricalColumn(i, groupLabelCatColor).custom('visible', show).width(70); // .asSet(',')
      builder.column(clustCol);
    } else if (col.metaInformation.listData) {
      // builder.column(LineUpJS.buildNumberColumn(i, [-10,10]).asArray().width(100));
      const columnDesc = deriveColumnDescriptions(data, { columns: [i] })[0];
      if (col.metaInformation.range) {
        columnDesc.domain = col.metaInformation.range;
      }

      if (col.metaInformation.colorMapping) {
        if (Array.isArray(col.metaInformation.colorMapping)) {
          columnDesc.colorMapping = {
            type: 'custom',
            entries: col.metaInformation.colorMapping.map((item, index) => {
              return { color: item, value: index / (col.metaInformation.colorMapping.length - 1) };
            }),
          };
        } else {
          columnDesc.colorMapping = col.metaInformation.colorMapping;
        }
        // column_desc["colorMapping"] = "interpolateBrBG";
      }
      builder.column(columnDesc);
    } else if (col.metaInformation.hide) {
      // don't show the column if e.g. it is only meta_data
    } else {
      builder.deriveColumns(i);
    }
  }

  // builder.deriveColumns([]);

  builder.column(LineUpJS.buildStringColumn('Annotations').editable().color(baseColor));
  builder.column(LineUpJS.buildStringColumn(UNIQUE_ID).width(50).color(baseColor)); // we need this to be able to filter by all indices; this ID corresponds to the mesh index

  builder.defaultRanking(true);
  // builder.deriveColors();
  builder.registerRenderer('mySmilesStructureRenderer', new MySmilesStructureRenderer());
  builder.registerRenderer('myLineChartRenderer', new MyLineChartRenderer());
  // builder.registerRenderer("myBarCellRenderer", new BarCellRenderer(true));
  builder.registerColumnType('mySmilesStructureColumn', StructureImageColumn);
  builder.registerColumnType('myLineChartColumn', TestColumn);
  builder.sidePanel(true, true); // collapse side panel by default
  builder.livePreviews({
    filter: false,
  });
  builder.dynamicHeight(myDynamicHeight);
  builder.animated(false);

  return builder;
}

export interface IStructureFilter extends IStringFilter {
  filter: string;
  valid: Set<string>;
}
export class StructureImageColumn extends StringColumn {
  protected structureFilter: IStructureFilter | null = null;

  filter(row: IDataRow): boolean {
    if (!this.isFiltered()) {
      return true;
    }
    return this.structureFilter!.valid.has(this.getLabel(row));
  }

  isFiltered(): boolean {
    return this.structureFilter != null && this.structureFilter.valid?.size > 0;
  }

  getFilter() {
    return this.structureFilter;
  }

  setFilter(filter: IStructureFilter | null) {
    if (isEqual(filter, this.structureFilter)) {
      return;
    }
    this.fire([StringColumn.EVENT_FILTER_CHANGED, Column.EVENT_DIRTY_VALUES, Column.EVENT_DIRTY], this.structureFilter, (this.structureFilter = filter));
  }
}

export class MyLineChartRenderer implements ICellRendererFactory {
  readonly title: string = 'Line Chart';

  canRender(col: TestColumn, mode: ERenderMode): boolean {
    // return col instanceof NumberColumn && (mode === ERenderMode.CELL);
    return mode === ERenderMode.CELL;
  }

  create(col: TestColumn): ICellRenderer {
    return {
      template: `<div class="svg-container">
        <svg class="svg-content" preserveAspectRatio="xMidYMid meet">
          <g>
            <path class="areaChart"></path>
            <path class="lineChart" fill="none" stroke="${baseColor}" stroke-width="0.02px"></path>
            <g>
              <line class="focus-line"></line>
              <line class="value-line-marker"></line>
            </g>
            <g>
              <text style="font-size:0.2px; opacity:0;" class="marker-text">-</text>
              <text style="font-size:0.2px;" class="focus-text"></text>
            </g>
            <rect class="hover-rect"></rect>
          </g>
        </svg></div>`,
      update: (n: HTMLImageElement, dataRow: IDataRow) => {
        if (renderMissingDOM(n, col, dataRow)) {
          return;
        }

        // get data
        const row = col.getMap(dataRow);
        const dataMeanList = row[0].value;
        const dataVarList = row[1].value;
        // const data_max = col.getMax();
        let dataMax =
          d3v5.max(dataMeanList, function (d) {
            return +d;
          }) + 1;
        // const data_min = col.getMin();
        let dataMin =
          d3v5.min(dataMeanList, function (d) {
            return +d;
          }) - 1;

        let measurementValue = null;
        let measurementStep = null;
        if (`${col.desc.label}_value` in dataRow.v && `${col.desc.label}_step` in dataRow.v) {
          measurementValue = dataRow.v[`${col.desc.label}_value`];
          measurementStep = dataRow.v[`${col.desc.label}_step`];

          dataMin = Math.min(dataMin, measurementValue - measurementValue * 0.1);
          dataMax = Math.max(dataMax, measurementValue + measurementValue * 0.1);
        }

        // this is the ratio that the chart should have
        const relWidth = WIDTH_HEIGHT_RATIO; // data_mean_list.length/4;
        const relHeight = 1;

        const div = d3v5.select(n);
        const svg = div.select('svg');
        svg.attr('viewBox', `0 0 ${relWidth} ${relHeight}`);

        // define x and y scales
        const x = d3v5.scaleTime().domain([0, dataMeanList.length]).range([0, relWidth]);

        const y = d3v5.scaleLinear().domain([dataMin, dataMax]).range([relHeight, 0]);

        // Show confidence interval
        svg
          .select('.areaChart')
          .datum(dataVarList)
          .attr('fill', '#c1c1c14d')
          .attr('stroke', 'none')
          .attr(
            'd',
            d3v5
              .area<number>()
              .x((d, i) => {
                return x(i);
              })
              .y0((d, i: number) => {
                return y(dataMeanList[i] - d);
              })
              .y1((d, i: number) => {
                return y(dataMeanList[i] + d);
              }),
          );

        // draw the line chart
        const path = svg.select('.lineChart');
        path.datum(dataMeanList).attr(
          'd',
          d3v5
            .line<number>()
            .x((d, i: number) => {
              return x(i);
            }) // i/data_list.length
            .y((d) => {
              return y(d);
            }), // 1-(d/data_max)
        );

        if (measurementValue != null && measurementStep != null) {
          // create the marker that marks an actual measurement

          svg
            .select('.value-line-marker')
            .style('fill', 'none')
            .attr('stroke', '#007dad')
            .attr('stroke-width', '0.5%')
            .attr('y1', '0')
            .attr('y2', relHeight)
            .attr('x1', x(measurementStep))
            .attr('x2', x(measurementStep))
            .style('opacity', 1);

          svg
            .select('.marker-text')
            .style('opacity', 1)
            .style('color', '#007dad')
            .attr('text-anchor', 'middle')
            .attr('alignment-baseline', 'middle')
            .attr('letter-spacing', '0px')
            .attr('x', x(measurementStep))
            .attr('y', y(measurementValue));
        }

        // add tooltips
        // https://www.d3-graph-gallery.com/graph/line_cursor.html
        // This allows to find the closest X index of the mouse:
        // var bisect = d3v5.bisector(function(d, i) { return i; }).left;

        // Create the line that travels along the x-axis of chart
        const focus = svg
          .select('.focus-line')
          .style('fill', 'none')
          .attr('stroke', 'black')
          .attr('stroke-width', '1%')
          .attr('y1', '0')
          .attr('y2', relHeight)
          .attr('x1', '0')
          .attr('x2', '0')
          .style('opacity', 0);

        // Create the text that travels along the curve of chart
        const focusText = svg
          .select('.focus-text')
          .style('opacity', 0)
          .attr('text-anchor', 'left')
          .attr('alignment-baseline', 'middle')
          .attr('letter-spacing', '0px');

        // What happens when the mouse move -> show the annotations at the right positions.
        function mouseover() {
          focus.style('opacity', 1);
          focusText.style('opacity', 1);
        }

        function mousemove() {
          // recover coordinate we need
          const x0 = d3v5.mouse(this)[0];
          // var y0 = d3v5.mouse(this)[1];
          // var i = bisect(data, x0, 1);
          // var x0 = x.invert(d3v5.mouse(this)[0]);
          let i = Math.round(x.invert(x0)); // x0*data_list.length

          i = Math.max(i, 0);
          i = Math.min(dataMeanList.length - 1, i);

          focus
            .attr('x1', x(i)) // i/data_list.length
            .attr('x2', x(i)); // i/data_list.length

          // // position the text in a way that it is always readable
          // if(x0 > rel_width/2){
          //     x0 = x0-rel_width/2;
          // }
          let measurementTxt = '';
          if (i === measurementStep) {
            measurementTxt = `<tspan x='0' dy='1.2em'>measured: ${Math.round(measurementValue * 100) / 100}</tspan>`;
          }
          focusText
            .html(
              `<tspan x='0' dy='1.2em'>step: ${i}</tspan><tspan x='0' dy='1.2em'>mean: ${
                Math.round(dataMeanList[i] * 100) / 100
              }</tspan><tspan x='0' dy='1.2em'>var: ${Math.round(dataVarList[i] * 100) / 100}</tspan>${measurementTxt}`,
            )
            .attr('x', 0) // x0
            .attr('y', 0); // y0
        }

        function mouseout() {
          focus.style('opacity', 0);
          focusText.style('opacity', 0);
        }

        // Create a rect on top of the svg area: this rectangle recovers mouse position
        svg
          .select('.hover-rect')
          .style('fill', 'none')
          .style('pointer-events', 'all')
          .attr('width', '100%')
          .attr('height', '100%')
          .on('mouseover', mouseover)
          .on('mousemove', mousemove)
          .on('mouseout', mouseout);
      },
    };
  }
}

export class MySmilesStructureRenderer implements ICellRendererFactory {
  readonly title: string = 'Compound Structure';

  // better with background image, because with img tag the user might drag the img when they actually want to select several rows
  readonly template = '<div style="background-size: contain; background-position: center; background-repeat: no-repeat;"></div>';

  canRender(col: StructureImageColumn, mode: ERenderMode): boolean {
    return col instanceof StructureImageColumn && (mode === ERenderMode.CELL || mode === ERenderMode.GROUP);
  }

  create(col: StructureImageColumn): ICellRenderer {
    return {
      template: this.template,
      update: (n: HTMLImageElement, d: IDataRow) => {
        // @ts-ignore
        const smiles = d.v[col.desc.column];
        ReactionCIMEBackendFromEnv.getStructureFromSmiles(smiles, false, null).then((x) => {
          if (x && x.length > 100) {
            // check if it is actually long enogh to be an img
            n.style.backgroundImage = `url('data:image/jpg;base64,${x}')`;
          } else {
            n.innerHTML = x;
          }
          n.title = `${mapSmilesToShortname(smiles)}: ${smiles}`;
          // n.alt = smiles;
        });
      },
    };
  }

  createGroup(col: StructureImageColumn, context: IRenderContext): IGroupCellRenderer {
    return {
      template: this.template,
      update: (n: HTMLImageElement, group: IOrderedGroup) => {
        const formData = new FormData();
        return context.tasks
          .groupRows(col, group, 'string', (rows) => {
            rows.every((row) => {
              const v = col.getLabel(row);
              formData.append('smiles_list', v);
              return true;
            });
          })
          .then(() => {
            ReactionCIMEBackendFromEnv.getMCSFromSmilesList(formData).then((x) => {
              n.style.backgroundImage = `url('data:image/jpg;base64,${x}')`;
              n.alt = formData.getAll('smiles_list').toString();
            });
          });
      },
    };
  }
}
