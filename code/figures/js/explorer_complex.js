// Define the data sources
var cplx_desc = cplx_desc_source.data;
var cplx_numeric = cplx_numeric_source.data;
var cplx_display = cplx_display_source.data;
var cplx_table = cplx_table_source.data
// Define the input widgets
var cplx_select = cplx_select_input.value;
var agg_method = agg_method_input.active;

// Define the aggregate options
var agg_options = ['min', 'max', 'median', 'mean']
var agg = agg_options[agg_method];

// Define the constant for the ecocyc address
ecocyc_link = "https://ecocyc.org/gene?orgid=ECOLI&id=" + cplx_select

// Find the empirical formula for the dataset and update display div and table
cplx_desc_field.text = '<b> EcoCyc Complex Identifier<b>: <a href=' + ecocyc_link + '><span style="font-family:Courier">' + cplx_select + '</span></a><br/><br/><b>Subunit Composition</b>'
var cplx_ind = getAllIndexes(cplx_desc['complex'], cplx_select);
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
console.log(minimum_subunits,cplx_table)
// Find the indices in the numeric dataset for the complex
var cplx_ind = getAllIndexes(cplx_numeric['complex'], cplx_select);


// Update the data. 
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

// Update the display source. 
cplx_display_source.change.emit()

// Custom Function Definitions
function getAllIndexes(arr, val) {
    var indices = [], i = -1;
    while ((i = arr.indexOf(val, i+1)) != -1){
        indices.push(i);
    }
    return indices;
}
