#%%
import numpy as np
import pandas as pd
import prot.viz
from bokeh.models import *
import bokeh.io
import bokeh.plotting
import scipy.stats
import bokeh.palettes
import bokeh.layouts
import tqdm
colors, palette = prot.viz.bokeh_theme()

# Load the data set and prune 
data = pd.read_csv('../../../data/compiled_annotated_complexes.csv', comment='#')
selector_df = pd.DataFrame([])
aggregate_df = []

# Define a unicode symbol lookup
unicode_digits = {0:'\u2080', 1:'\u2081', 2:'\u2082', 3:'\u2083', 4:'\u2084',
                  5:'\u2085', 6:'\u2086', 7:'\u2087', 8:'\u2088', 9:'\u2089'}
for i in range(10, 1000):
    ints = list(str(i))
    entry = ''
    for k in ints:
        entry += f'{unicode_digits[int(k)]}'
    unicode_digits[i] = entry

#%%
# Define the markers for the data set
markers = {'li_2014':'square', 'schmidt_2016':'circle', 
            'valgepea_2013':'diamond', 'peebo_2015':'triangle'} 

for g, d in tqdm.tqdm(data.groupby(['complex']), desc='Iterating through complexes'):
    # Get the gene constituent and compute the stoichiometric formula
    formula = ''
    _d = d.drop_duplicates(subset=['gene_name', 'n_subunits'])
    genes = _d['gene_name'].values
    n_units = _d['n_subunits'].values
    for gene, n_units in zip(genes, n_units):
        n_units = str(int(n_units))
        if '-' in gene:
            gene_name = f'{gene[0].upper()}{gene[1:-3]}{gene[-3].upper()}{gene[-2:]}'
        else:
            gene_name = f'{gene[0].upper()}{gene[1:-1]}{gene[-1].upper()}'
        formula +=  f'[{gene_name}]{unicode_digits[int(n_units)]}'

    # Compute min, max, mean, and median number of subunits and the
    # corresponding spearman R
    pooled_stats = d.groupby(['growth_rate_hr'])['n_units'].agg(('mean', 'median', 'min', 'max')).reset_index()
    pooled_stats.dropna(inplace=True)
    correlation = [scipy.stats.spearmanr(pooled_stats['growth_rate_hr'], pooled_stats[stat])[0] for stat in ['mean', 'median', 'min', 'max']]
    for _g, _d in d.groupby(['dataset', 'dataset_name', 'condition']): 
        stats = _d.groupby(['growth_rate_hr'])['n_units'].agg(
                    ('mean', 'median', 'min', 'max')).reset_index()
        # Fix gross naming in the annotated function. 
        function = d['complex_annotation'].values[0]
        function = function.replace('<i>C</i>', 'C')
        function = function.replace('<i>N</i>', 'N')
        function = function.replace('<i>O</i>', 'O')
        function = function.replace('<sup>', '')
        function = function.replace('<sup>+</sup>', '\u207A')
        function = function.replace('<sup>-</sup>', '\u207B')
        function = function.replace('<sup>2-</sup>', '\u00B2\u207B')
        function = function.replace('<sup>2+</sup>', '\u00B2\u207A')
        function = function.replace('<sup>3+</sup>', '\u00B3\u207A')
        function = function.replace('<sup>4+</sup>', '\u2084\u207A')
        function = function.replace('</sup>', '')
        function = function.replace('<sub>', '')
        function = function.replace('</sub>', '')
        function = function.replace('&alpha;', 'α')
        function = function.replace('&beta;', 'β')
        function = function.replace('&gamma;', 'γ')
        function = function.replace('&mdash;', '-')


        stats['marker'] = markers[_g[0]]
        stats['dataset'] = _g[1]
        stats['growth_condition'] = _g[-1]
        stats['complex'] = g
        stats['function'] = function

        # Append the aggregate dataframe
        aggregate_df.append(stats)

    # Append the name information to the selector dataframe
    selector_df = selector_df.append({'cplx':g, 
                      'formula':formula,
                      'function':function,
                      'cog_category':d['cog_category'].values[0],
                      'cog_class':d['cog_class'].values[0],
                      'r_mean':correlation[0],
                      'r_median':correlation[1],
                      'r_min':correlation[2],
                      'r_max':correlation[3]},
                      ignore_index=True)

aggregate_data = pd.concat(aggregate_df, sort=False)
# %%
aggregate_data.dropna(inplace=True)
bokeh.io.output_file('../../../figures/complex_explorer.html')

# Define the selection source and widget
selection_options = [(b, a) for a, b in zip(selector_df['function'].unique(), selector_df['cplx'].unique())]

# Define the display sources. 
displayed = ColumnDataSource({'x':[], 'y':[], 'c':[], 'fn':[], 'm':[],
                               'ds':[], 'gc':[]})
source = ColumnDataSource(aggregate_data)
color_source = np.array(bokeh.palettes.Category20_20)
np.random.shuffle(color_source)

# Define the interaction widgets
selection = MultiChoice(options=selection_options, solid=False, title='displayed complexes') 
statistic = RadioButtonGroup(labels=['mean', 'median', 'minimum', 'maximum'], active=0)

# Define the divs
desc = Div(text="Description of displayed complexes will be reported here")

# Set up the plotting axis
canvas = bokeh.plotting.figure(width=600, height=400, x_axis_label='growth rate [hr^-1]',
                               y_axis_label='possible complexes per cell',
                               x_range=[0, 2])

canvas.scatter(x='x', y='y', fill_color='c', marker='m', size=10, source=displayed)
# 
cb = """
var displayedData = DisplaySource.data;
var sourceData = DataSource.data;
var complexes = SelectionWidget.value;
var statistic = StatisticWidget.active;

// Update the plot labels
var statisticLabels = ['average', 'median', 'minimum', 'maximum'];
var statisticKeys = ['mean', 'median', 'min', 'max']
// iterate through the complexes
var x = []; // Growth rate
var y = []; // Complex count
var c = []; // Complex color
var m = []; // Marker style
var ds = []; // Dataset name
var fn = []; // Complex function
var gc = []; // Growth condition
for (var i = 0; i < complexes.length; i++) {
    var selCplx = complexes[i];
    for (var j = 0; j < sourceData.complex.length; j++) {
        if (sourceData.complex[j] === selCplx) { 
            x.push(sourceData['growth_rate_hr'][j]);
            y.push(sourceData[statisticKeys[statistic]][j]);
            c.push(colors[i]);
            m.push(sourceData['marker'][j]);
            ds.push(sourceData['dataset'][j]);
            fn.push(sourceData['function'][j]);
            gc.push(sourceData['growth_condition'][j]);
            
        }
    }

}
console.log(displayedData)
displayedData['x'] = x;
displayedData['y'] = y;
displayedData['c'] = c;
displayedData['m'] = m;
displayedData['fn'] = fn;
displayedData['ds'] = ds;
displayedData['gc'] = gc;
DisplaySource.change.emit()
"""

js = CustomJS(args={'DisplaySource':displayed, 'DataSource':source,
                    'StatisticWidget':statistic, 'SelectionWidget':selection,
                    'canvas':canvas, 'colors':color_source},
             code=cb)

selection.js_on_change('value', js)
statistic.js_on_change('active', js)
col = bokeh.layouts.column(statistic, canvas, selection)
row = bokeh.layouts.row(col, desc)
lay = row
bokeh.io.save(lay)

# %%
