// Define the data sources
var cplx_desc = cplx_desc_source.data;
var cplx_numeric = cplx_numeric_source.data;
var cplx_display = cplx_display_source.data;
var cplx_table = cplx_table_source.data

// Define the input values
var agg_method = agg_method_input.active;
var agg_options = ['min', 'max', 'median', 'mean']
var agg = agg_options[agg_method];

// Define the input widgets
var cplx_annotation = cplx_select_input.value;

// Define the constant for the ecocyc address
var cplx_ind = getAllIndexes(cplx_desc['complex_annotation'], cplx_annotation);
var cplx_select = cplx_desc['complex'][cplx_ind[0]];


// Update the display table
ecocyc_link = "https://ecocyc.org/gene?orgid=ECOLI&id=" + cplx_select
cplx_desc_field.text = '<b> EcoCyc Complex Identifier<b>: <a href=' + ecocyc_link + '><span style="font-family:Courier">' + cplx_select + '</span></a><br/><br/><b>Subunit Composition</b>'

for (var i = 0; i < cplx_ind.length; i++) {
    if (i == 0) {
        cplx_table['protein'].length = 0;
        cplx_table['subunits'].length = 0;
        cplx_table['func'].length = 0;
        cplx_table['relative_subunits'].length = 0;
    }
    cplx_table['protein'].push(cplx_desc['protein'][cplx_ind[i]]);
    cplx_table['subunits'].push(cplx_desc['subunits'][cplx_ind[i]]);
    cplx_table['func'].push(cplx_desc['func'][cplx_ind[i]]);
}
var minimum_subunits = Math.min(...cplx_table['subunits']);
for (var i = 0; i < cplx_ind.length; i++) {
    cplx_table['relative_subunits'].push(cplx_table['subunits'][i] / minimum_subunits);
    cplx_table['observed'].push('hover to view');
} 
cplx_table_source.change.emit();

// Update the data. 
cplx_ind = getAllIndexes(cplx_numeric['complex_annotation'], cplx_annotation);
for (var i = 0; i < cplx_ind.length; i++) {
    if (i == 0) { 
        cplx_display['x'].length = 0;
        cplx_display['y'].length = 0;
        cplx_display['c'].length = 0;
        cplx_display['condition'].length = 0;
        cplx_display['l'].length = 0;
    }
    cplx_display['x'].push(cplx_numeric['growth_rate_hr'][cplx_ind[i]]);
    cplx_display['y'].push(cplx_numeric[agg][cplx_ind[i]]);
    cplx_display['c'].push(cplx_numeric['color'][cplx_ind[i]]);
    cplx_display['l'].push(cplx_numeric['dataset_name'][cplx_ind[i]]);
    cplx_display['condition'].push(cplx_numeric['condition'][cplx_ind[i]]);
    cplx_display['cond'].push(cplx_numeric['cond'][cplx_ind[i]]);
}
cplx_display_source.change.emit()

// Custom Function Definitions
function getAllIndexes(arr, val) {
    var indices = [], i = -1;
    while ((i = arr.indexOf(val, i+1)) != -1){
        indices.push(i);
    }
    return indices;
}
