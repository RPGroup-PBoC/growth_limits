var display = complex_display.data;
var source = complexes.data;
var complex = complex_selector.value;
var agg = agg_fn.active;
var agg_opts = ['min', 'max', 'median', 'mean']
var key = agg_opts[agg];
var complex_inds = getAllIndexes(source['complex_annotation'], complex)
for (var i = 0; i < complex_inds.length; i++) {
    if (i == 0) { 
        display['x'].length = 0;
        display['y'].length = 0;
        display['c'].length = 0;
        display['l'].length = 0;
        display['cond'].length = 0;
    }
    display['x'].push(source['growth_rate_hr'][complex_inds[i]]);
    display['y'].push(source[key][complex_inds[i]]);
    display['l'].push(source['dataset_name'][complex_inds[i]]);
    display['c'].push(source['color'][complex_inds[i]]);
    display['cond'].push(source['condition'][complex_inds[i]])
}
console.log(display)
complex_display.change.emit();

// Update the display table
var ecocyc_inds = getAllIndexes(desc.data['complex_annotation'], complex)
var ecocyc_cplx = desc.data['complex'][ecocyc_inds[0]]
ecocyc_link = "https://ecocyc.org/gene?orgid=ECOLI&id=" + ecocyc_cplx
complex_desc.text = '<b> EcoCyc Complex Identifier<b>: <a target="_blank" href=' + ecocyc_link + '><span style="font-family:Courier">' + ecocyc_cplx + '</span></a><br/>'


// Custom Function Definitions
function getAllIndexes(arr, val) {
    var indices = [], i = -1;
    while ((i = arr.indexOf(val, i+1)) != -1){
        indices.push(i);
    }
    return indices;
}