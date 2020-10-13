var display = protein_display.data;
var source = proteins.data;
var protein = protein_selector.value;
var prot_inds = getAllIndexes(source['protein'], protein);
var product = source['gene_product'][prot_inds[0]];
for (var i = 0; i < prot_inds.length; i++) {
    if (i == 0) { 
        display['x'].length = 0;
        display['y'].length = 0;
        display['c'].length = 0;
        display['l'].length = 0;
        display['cond'].length = 0;
    }
    display['x'].push(source['growth_rate_hr'][prot_inds[i]]);
    display['y'].push(source['abundance'][prot_inds[i]]);
    display['l'].push(source['dataset_name'][prot_inds[i]]);
    display['c'].push(source['color'][prot_inds[i]]);
    display['cond'].push(source['condition'][prot_inds[i]])
}

protein_display.change.emit();
protein_desc.text = "<b>EcoCyc Gene Product Annotation</b>:<br/>" + product

// Custom Function Definitions
function getAllIndexes(arr, val) {
    var indices = [], i = -1;
    while ((i = arr.indexOf(val, i+1)) != -1){
        indices.push(i);
    }
    return indices;
};