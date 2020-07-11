#%%
import numpy as np
import bokeh.io
import bokeh.plotting
from bokeh.models import * 
import bokeh.layouts
import prot.estimate
import prot.viz 
colors, palette = prot.viz.bokeh_theme()
bokeh.io.output_file('model_explorer.html')


# Define the data sources
ribosomes = np.logspace(3, 6, 200)
elongation_rate = [0]
data_source = ColumnDataSource({'R': ribosomes, 
                                'R_fa': ribosomes,
                                'elongation_rate':np.zeros_like(ribosomes), 
                                'elongation_rate_fa':np.zeros_like(ribosomes), 
                                'lambda': np.zeros_like(ribosomes),
                                'leg': ['100% ribosomes active'],
                                'lambda_fa': np.zeros_like(ribosomes),
                                'supply': [1E5],
                                'consumed': elongation_rate,
                                'consumed_fa': elongation_rate,
                                'accum': [0],
                                'accum_fa':[0]})
                                


# Define the sliders. 
kd_slider = Slider(value=5, start=0.01, end=100, step=0.01, title='Kd [mM]')
r_aa = Slider(value=5, start=3, end=8, step=0.1, title='log\u2081\u2080 (amino acid supply rate [AA / (s•µm\u00b3)])')
fa_slider = Slider(value=1, start=0.001, end=1, step=0.001, title='active ribosome fraction')


# Set up the figure canvases. 
elong_rate_plot = bokeh.plotting.figure(width=500, height=500, 
                                      x_axis_label='ribosomes per µm\u00b3',
                                      x_axis_type='log',
                                      y_axis_label='elongation rate [AA / sec]',
                                      x_range=[1E3, 5E5],
                                      y_range = [0, 18])

growth_rate_plot = bokeh.plotting.figure(width=500, height=500, 
                                      x_axis_label='ribosomes per µm\u00b3',
                                      x_axis_type='log',
                                      y_axis_label='growth rate [hr^-1]',
                                      x_range=[1E3, 5E5],
                                      y_range = [0, 2])

conc_plot = bokeh.plotting.figure(width=500, height=500, 
                                      x_axis_label='ribosomes per µm\u00b3',
                                      x_axis_type='log',
                                      y_axis_type='log',
                                      y_axis_label='AA / (s • µm\u00b3)',
                                      x_range=[1E3, 5E5],
                                      y_range=[1E0, 1E9])





# Elongation Plot
elong_rate_plot.line(x='R', y='elongation_rate_fa', line_color=colors['light_blue'],
                     source=data_source, line_width=3, legend_field='leg')

elong_rate_plot.line(x='R', y='elongation_rate', line_color=colors['blue'],
                     source=data_source, line_width=3, legend_label='no ribosome regulation')

# Growth Rate Plot
growth_rate_plot.line(x='R', y='lambda_fa',  line_width=3,
                     source=data_source, line_color=colors['light_blue'],
                      legend_field='leg')

growth_rate_plot.line(x='R', y='lambda',  line_width=3, legend_label='no ribosome regulation',
                     source=data_source, line_color=colors['blue'])

# Concentration Plot
conc_plot.line(x='R', y='supply', line_width=3, legend_label='[AA] generated',
               line_color='black', source=data_source)
conc_plot.line(x='R', y='consumed_fa', line_width=3, legend_label='[AA] consumed (with regulation)',
              line_color=colors['light_red'], source=data_source)
conc_plot.line(x='R', y='consumed', line_width=3, legend_label='[AA] consumed (no regulation)',
              line_color=colors['red'], source=data_source)
conc_plot.line(x='R', y='accum_fa', line_width=3, legend_label='[AA]_eff (with regulation)',
              line_color=colors['light_purple'], source=data_source)
conc_plot.line(x='R', y='accum', line_width=3, legend_label='[AA]_eff (no regulation)',
              line_color=colors['purple'], source=data_source)


# Load the javascript. 
with open('model_explorer.js') as f:
    js = f.read()
cb = CustomJS(args={'source':data_source, 'kd_slider':kd_slider,  'raa_slider':r_aa, 'fa_slider': fa_slider},
              code=js)

# Assign the callbacks
kd_slider.js_on_change('value', cb)
r_aa.js_on_change('value', cb)
fa_slider.js_on_change('value', cb)
kd_slider.callback = cb
r_aa.callback = cb
fa_slider.callback = cb

# Define the layout
controls = bokeh.layouts.row(kd_slider, r_aa, fa_slider)
plots = bokeh.layouts.row(elong_rate_plot, growth_rate_plot, conc_plot)
lay = bokeh.layouts.column(controls, plots)
bokeh.io.save(lay)


# %%


# %%
