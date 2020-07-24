var prot_source = prot_numeric_source.data;
var prot_display = prot_display_source.data;
var prot_select = prot_select_input.value;

// Find the indices of the selected protein
var prot_inds = getAllIndexes(prot_source['gene_name'], prot_select)
console.log(prot_display)
for (var i = 0; i < prot_inds.length; i++) {
    if (i == 0) { 
        prot_display['x'].length = 0;
        prot_display['y'].length = 0;
        prot_display['c'].length = 0;
        prot_display['l'].length = 0;
        prot_display['condition'].length = 0;
    }

    prot_display['x'].push(prot_source['growth_rate_hr'][prot_inds[i]]);
    prot_display['y'].push(prot_source['tot_per_cell'][prot_inds[i]]);
    prot_display['c'].push(prot_source['color'][prot_inds[i]]);
    prot_display['l'].push(prot_source['dataset_name'][prot_inds[i]]);
    prot_display['condition'].push(prot_source['condition'][prot_inds[i]]);
}

prot_display_source.change.emit();


// Get the indices where the condition is met. 
function getAllIndexes(arr, val) {
    var indices = [], i = -1;
    while ((i = arr.indexOf(val, i+1)) != -1){
    indices.push(i);
    }
    return indices;
}