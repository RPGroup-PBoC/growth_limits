// Identify the provided selected class
var class_select = class_select_input.value;

// Given class selection, restrict selectable complexes to those in annotated class
var cplx_indices = getAllIndexes(cplx_desc_source.data['cog'], class_select);
var prot_indices = getAllIndexes(prot_desc_source.data['cog'], class_select);
console.log(prot_indices)
var cplx_menu = [];
for (var i = 0; i < cplx_indices.length; i++ ) {
    var annotation = cplx_desc_source.data['complex_annotation'][cplx_indices[i]];
    // Determine if the annotation is already registered in the dropdown. 
    var menu_inds = getAllIndexes(cplx_menu, annotation);
    if (menu_inds.length == 0) {
        cplx_menu.push(annotation);
    }
    cplx_menu.sort()
}
var prot_menu = [];
for (var i = 0; i< prot_indices.length; i++) {
    var prot = prot_desc_source.data['gene_name'][prot_indices[i]];
    prot_menu.push(prot);
}
prot_menu.sort()
// Update the menu on the selector
complex_select_input.options =  cplx_menu;
complex_select_input.value = cplx_menu[0];
complex_select_input.change.emit();
prot_select_input.options = prot_menu;
prot_select_input.value = prot_menu[0];
prot_select_input.change.emit();

// Custom Function Definitions
function getAllIndexes(arr, val) {
    var indices = [], i = -1;
    while ((i = arr.indexOf(val, i+1)) != -1){
        indices.push(i);
    }
    return indices;
}
