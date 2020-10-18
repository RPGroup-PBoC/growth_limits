// Define variables
var kd = kd_slider.value * 0.001;
var r_aa = Math.pow(10, raa_slider.value);
var fa = fa_slider.value;

// Define constants. 
const rt_max = 17.1; // in AA  / sec 
const avo = 6.022E23; // Avogadro's number
const V = 1E-15 ; // VOlume in liters

// Define data sources
var data = source.data;
var R = data['R'];
var Ra = [];
var label = Number(fa * 100).toFixed(1)  + '% active';
var leg = [];
var supply = [];
var consumed = [];
var consumed_fa = [];
var accum = [];
var accum_fa = [];

for (var i = 0; i < R.length; i++) {
    Ra[i] = fa * R[i];
    leg[i] = label;
}


var NAA = compute_NAA(R);
var vol = compute_vol(R);

// Given updated values, compute new values. 
var r_t = compute_elong_rate(r_aa, R, kd, 1, V); 
var r_t_fa = compute_elong_rate(r_aa, R, kd, fa, V);
var lambda = compute_growth_rate(NAA, R, 1, vol);
var lambda_fa = compute_growth_rate(NAA, R, fa, vol);
data['R_fa'] = Ra;
data['lambda'] = lambda;
data['lambda_fa'] = lambda_fa;
data['elongation_rate'] = r_t;
data['elongation_rate_fa'] = r_t_fa;
data['leg'] = leg;

// Compute the consumption and production levels. 
var R_cell = [];
for (var i = 0; i <  R.length; i++) {
    supply[i] = r_aa;
    consumed[i] = r_t[i] * R[i];
    consumed_fa[i] = r_t[i] * R[i] * fa;
    accum[i] = r_aa - r_t[i] * R[i];
    accum_fa[i] = r_aa - r_t[i] * R[i] * fa;
    R_cell[i] = vol[i] * R[i];
}
data['R_cell'] = R_cell;
data['supply'] = supply;
data['consumed'] = consumed;
data['consumed_fa'] = consumed_fa;
data['accum'] = accum;
data['accum_fa'] = accum_fa;

// Push the changes to the DOM.
source.change.emit();

// -----------------------------------------------------------------------------
// Function Definition
// -----------------------------------------------------------------------------

function compute_elong_rate(r_aa, R, Kd, fa) {
    var r_t = [];
    for (var i = 0; i < R.length; i++) {
        var a = -R[i] * fa;
        var b = r_aa + rt_max * R[i] * fa + Kd * V * avo;
        var c = -rt_max * r_aa;
        var numer = -b + Math.sqrt(Math.pow(b, 2) - 4 * a * c);
        r_t[i] = numer / (2 * a);
    }
    return r_t
}

function compute_elong_rate_volume(r_aa, R, Kd, fa, volume) {
    var r_t = [];
    for (var i = 0; i < R.length; i++) {
        var a = -R[i] * fa;
        var b = r_aa + rt_max * R[i] * fa + Kd * volume[i] * 1E-15 * avo;
        var c = -rt_max * r_aa;
        var numer = -b + Math.sqrt(Math.pow(b, 2) - 4 * a * c);
        r_t[i] = numer / (2 * a);
    }
    return r_t
}

function compute_growth_rate(Naa, R, fa, volume) {
    var lam = [];
    var v_rt = compute_elong_rate_volume(r_aa, R, kd, fa, volume);
    for (var i = 0; i < R.length; i++ ) { 
        lam.push(3600 * v_rt[i] * R[i] * fa / Naa[i])
    }
    return lam
}

function compute_NAA(R) { 
    var NAA = []
    for (var i = 0 ; i < R.length; i++) { 
    NAA.push(Math.exp(16.1) * Math.pow(R[i], 0.5));
    }
    return NAA
}

function compute_vol(R) { 
    var vol = []
    for (var i = 0; i < R.length; i++) {
    vol.push(Math.exp(-5.8) * Math.pow(R[i], 0.6));
    }
    return vol
}