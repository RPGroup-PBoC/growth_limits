// Get the hover index
const hover_idx = cb_data.index.indices 
var cplx_source = cplx_source_data.data;
var cplx_select = cplx_select_input.value;
var cplx_table = cplx_table_source.data;
var cplx_display = cplx_display_source.data;
var condition = cplx_display['cond'][hover_idx];
var dataset = cplx_display['l'][hover_idx];
var hr = cplx_display['x'][hover_idx];

// Find the indices of the complex 
var _cplx_inds = getAllIndexes(cplx_desc_source.data['complex_annotation'], cplx_select)
var cplx = cplx_desc_source.data['complex'][_cplx_inds[0]]
var cplx_inds = getAllIndexes(cplx_source['complex'], cplx);
var approved_inds = [];
for (var i = 0; i < cplx_inds.length; i++) {
    if ((cplx_source['dataset_name'][cplx_inds[i]]==dataset) & (cplx_source['growth_rate_hr'][cplx_inds[i]]==hr)) {
        approved_inds.push(cplx_inds[i]);
    }
}
if (approved_inds.length == 0) { 
    cplx_table['observed'].length = 0;
    for (var i = 0; i < cplx_table['subunits'].length; i++) { 
        cplx_table['observed'].push('hover to view');
    }
        
 cplx_table_source.change.emit();
}

else { 
    // Update the table
    var totals = []
    for (var i = 0; i < cplx_table['subunits'].length; i++) {
        // Find the indices of the gene
        for (var j = 0; j < approved_inds.length; j++) {
            if (cplx_source['gene_name'][approved_inds[j]].toLowerCase()==cplx_table['protein'][i].toLowerCase()) {
                totals.push(cplx_source['tot_per_cell'][approved_inds[j]]);
            }
        }
    }
    var min_ind = cplx_table['relative_subunits'].indexOf(1); 
    cplx_table['observed'].length = 0;
    for (var i = 0; i < totals.length; i++) {
        var val = totals[i] / totals[min_ind]
        var val = val.toFixed(2);
        cplx_table['observed'].push(val);
      } 
    cplx_table_source.change.emit();
}
// Get the indices where the condition is met. 
function getAllIndexes(arr, val) {
    var indices = [], i = -1;
    while ((i = arr.indexOf(val, i+1)) != -1){
    indices.push(i);
    }
    return indices;
}


