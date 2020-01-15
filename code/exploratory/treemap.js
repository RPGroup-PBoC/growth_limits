d3 = require()

// Set the dimensions of the plot
var margin = {top:10, bottom:10, right: 10, left:10};
var width = 445 - margin.left - margin.right;
var height = 445 - margin.top - margin.bottom;


// Add the SVG object to teh page. 
var svg = d3.select("#viz")
    .append("svg")
    .attr("width", width)
    .attr("height", height)
    .append("g")
    .attr("transform", "translate(" + margin.left + ',' + margin.top + ")");

// Read in the data
// d3.csv('../../')
