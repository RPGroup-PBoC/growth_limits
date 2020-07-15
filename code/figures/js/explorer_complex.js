// Define the data sources
var cplx_desc = cplx_desc_source.data;
var cplx_numeric = cplx_numeric_source.data;

// Define the input widgets
var cplx_select = cplx_select_input.value;
var scaling = scaling_input.value;
var agg_method = agg_method_input.value;

console.log(cplx_select, scaling, agg_method)

// Define the constant for the ecocyc address
ecocyc_link = "https://ecocyc.org/gene?orgid=ECOLI&id=" + cplx_select[1]

// Find the empirical formula for the dataset
// for (var i = 0; i< cplx_desc.length; i++) {
var cplx_ind = cplx_desc['complex'].indexOf(cplx_select);
var formula = cplx_desc['formula'][cplx_ind];

// Update the outputdiv
cplx_desc_field.text = formula + '<br/> <a href=${ecocyc_link}${cplx_select}>EcoCyc Entry</a>';


