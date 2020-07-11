// Define variables
var kd = kd_slider.value;
var r_aa = Math.pow(10, raa_slider.value);
var fa = fa_slider.value;


// Define constants. 
var rt_max = 17.1; // in AA  / sec 
var avo = 6.022E5; // Avogadro's number for conversion given milimolar concs at fL vol
var rho = 1.1 // cell density in pg / fl
var theta_prot =  0.15; // Fraction of dry mass that is protein
var m_aa = 110 / 6.022E11; // Mass of average amino acid in pg
var vol = 1 // in µm^3
var Lr = 7500 // AA per ribosome

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


console.log(Ra)
// Given updated values, compute new values. 
var r_t = compute_elong_rate(R, r_aa, kd, vol, 1); 
var r_t_fa = compute_elong_rate(R, r_aa, kd, vol, fa);
var lambda = compute_growth_rate(r_t, R, vol, 1);
var lambda_fa = compute_growth_rate(r_t, R, vol, fa);
data['R_fa'] = Ra;
data['lambda'] = lambda;
data['lambda_fa'] = lambda_fa;
data['elongation_rate'] = r_t;
data['elongation_rate_fa'] = r_t_fa;
data['leg'] = leg;

// Compute the consumption and production levels. 
for (var i = 0; i <  R.length; i++) {
    supply[i] = r_aa;
    consumed[i] = r_t[i] * R[i];
    consumed_fa[i] = r_t[i] * R[i] * fa;
    accum[i] = r_aa - r_t[i] * R[i];
    accum_fa[i] = r_aa - r_t[i] * R[i] * fa;
}
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

function compute_NAA(vol) { 
    // Computes the number of amino acids per cell based on a volume
    return (rho * vol * theta_prot) / (m_aa);
}

function compute_elong_rate(R, r_aa, kd, vol, pact) { 
    // Computes the quadratic solution for the elongation rate. 

    // Compute the various components. 
    var NAA = compute_NAA(vol);
    var r_t = [];
    for (var i = 0; i < R.length; i++) {
        var theta = (kd * vol * avo) / (Math.log(2) * (NAA + Lr * R[i]));
        var root = Math.sqrt(4 * theta * rt_max * r_aa * R[i] * pact + Math.pow(rt_max * R[i] * pact, 2) - 2 * rt_max * r_aa * R[i] * pact + Math.pow(r_aa, 2));
        r_t[i] = (rt_max * R[i] * pact + r_aa - root) / (2 * R[i] * pact * (1 - theta)); 
    }
    return r_t
}

function compute_growth_rate(r_t, R, vol, pact) {
    // Given  an elongation rate and a unit volume, compute the growth rate
    var NAA = compute_NAA(vol);
    var lambda = [];
    for (var i = 0; i < R.length; i++) {
        lambda[i] = 3600 * (r_t[i] * R[i] * pact)  / (NAA + R[i] * Lr); 
    }
    return lambda;
}
