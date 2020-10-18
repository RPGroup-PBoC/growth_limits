// Identify the provided selected class
var cog = COG_selector.value;

// Given class selection, restrict selectable complexes to those in annotated class
var cog_indices = getAllIndexes(desc.data['class'], cog);
var cplx_menu = [];
for (var i=0; i < cog_indices.length; i++) { 
    cplx_menu.push(desc.data['complex_annotation'][cog_indices[i]]);
}

cplx_menu.sort()
complex_selector.options =  cplx_menu;
complex_selector.value = cplx_menu[0];
complex_selector.change.emit();

// Custom Function Definitions
function getAllIndexes(arr, val) {
    var indices = [], i = -1;
    while ((i = arr.indexOf(val, i+1)) != -1){
        indices.push(i);
    }
    return indices;
}
