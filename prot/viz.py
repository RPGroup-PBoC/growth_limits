import pandas as pd
import squarify as sq
from bokeh.themes import Theme
import matplotlib.pyplot as plt
import seaborn as sns
import altair as alt
import bokeh.io

def assign_rect_bounds(df, key, width=500, height=500, text_pad=0,
                       copy_df=False):
    """
    Uses squarify to assign the bounds for a tree map. 
    """
    # Make a copy of the dataframe
    if copy_df:
        df = df.copy()

    # Compute the values and rects. 

    # Drop the quantities that are equal to zero in hte key dimension
    df = df[df[key] > 0]
    values = sq.normalize_sizes(df[key], width, height)
    rects = sq.squarify(values, 0, 0, width, height)

    # Assign the information to the data frame and return
    df['left'] = [r['x'] for r in rects] 
    df['right'] = [r['x'] + r['dx'] for r in rects]
    df['bottom'] = [r['y'] for r in rects]
    df['top'] = [r['y'] + r['dy'] for r in rects]

    # Determine if text padding should be applied
    df['text_left'] = df['left'] + text_pad
    df['text_bottom'] = df['bottom'] + text_pad
    df['text_right'] = df['right'] - text_pad
    df['text_top'] = df['top'] - text_pad

    return df

def altair_theme():
    colors = {'green': '#7AA974', 'light_green': '#BFD598',
              'pale_green': '#DCECCB', 'yellow': '#EAC264',
              'light_yellow': '#F3DAA9', 'pale_yellow': '#FFEDCE',
              'blue': '#738FC1', 'light_blue': '#A9BFE3',
              'pale_blue': '#C9D7EE', 'red': '#D56C55', 'light_red': '#E8B19D',
              'pale_red': '#F1D4C9', 'purple': '#AB85AC',
              'light_purple': '#D4C2D9', 'dark_green':'#7E9D90', 'dark_brown':'#905426'}
    palette = [colors['red'], colors['blue'], colors['green'], 
               colors['purple'], colors['dark_green'], colors['dark_brown'],
               colors['yellow']]
    
    def _theme():
        return {
            'config': {
                'background': 'white',
                    'group': { 
                    'fill': '#E3DCD0'
                    },
                'view': {
                    'strokeWidth': 0,
                    'height': 300,
                    'width': 400,
                    'fill': '#E3DCD0'
                    },
                'mark': {
                    'strokeWidth': 0.5,
                    'stroke': 'black'
                },
                'axis': {
                    'domainColor': None,
                    'labelFont': 'Lucida Sans',
                    'titleFont': 'Lucida Sans',
                    'titleFontWeight': 400,
                    'grid': False,
                    'ticks': True,
                    'tickColor': 'white',
                    'tickOffset': 8,
                    'tickWidth': 1.5
                },
                'range': {
                    'category': palette
                },
                'legend': {
                    'labelFont': 'Lucida Sans',
                    'titleFont': 'Lucida Sans',
                    'titleFontWeight': 400
                },
                'title' : { 
                    'font': 'Lucida Sans',
                    'fontWeight': 400,
                    'anchor': 'middle'

                }
                  }
                }

    alt.themes.register('pboc', _theme)# enable the newly registered theme
    alt.themes.enable('pboc')



def bokeh_theme():
    """A custom bokeh theme to match PBoC 2e colors"""
    theme_json = {'attrs':
            {'Figure': {
                'background_fill_color': '#E3DCD0',
                'outline_line_color': '#FFFFFF',
            },
            'Axis': {
            'axis_line_color': "white",
            'major_tick_in': 7,
            'major_tick_line_width': 2.5,
            'major_tick_line_color': "white",
            'minor_tick_line_color': "white",
            'axis_label_text_font': 'Helvetica',
            'axis_label_text_font_style': 'normal'
            },
            'Grid': {
                'grid_line_color': None,
            },
            'Legend': {
                'background_fill_color': '#E3DCD0',
                'border_line_color': '#FFFFFF',
                'border_line_width': 1.5,
                'background_fill_alpha': 0.5
            },
            'Text': {
                'text_font_style': 'normal',
               'text_font': 'Helvetica'
            },
            'Title': {
                'background_fill_color': '#FFEDC0',
                'text_font_style': 'normal',
                'align': 'center',
                'text_font': 'Helvetica',
                'offset': 2,
            }}}

    theme = Theme(json=theme_json)
    bokeh.io.curdoc().theme = theme

    # Define the colors
    colors = {'green': '#7AA974', 'light_green': '#BFD598',
              'pale_green': '#DCECCB', 'yellow': '#EAC264',
              'light_yellow': '#F3DAA9', 'pale_yellow': '#FFEDCE',
              'blue': '#738FC1', 'light_blue': '#A9BFE3',
              'pale_blue': '#C9D7EE', 'red': '#D56C55', 'light_red': '#E8B19D',
              'pale_red': '#F1D4C9', 'purple': '#AB85AC',
              'light_purple': '#D4C2D9', 'dark_green':'#7E9D90', 'dark_brown':'#905426'}
    palette = [v for k, v in colors.items() if 'pale' not in k]
    return [colors, palette]


def plotting_style(grid=False):
    """
    Sets the style to the publication style
    """
    rc = {'axes.facecolor': '#E3DCD0',
          'font.family': 'Lucida Sans Unicode',
          'grid.linestyle': '-',
          'grid.linewidth': 0.5,
          'grid.alpha': 0.75,
          'grid.color': '#ffffff',
          'axes.grid': grid,
          'ytick.direction': 'in',
          'xtick.direction': 'in',
          'xtick.gridOn': grid,
          'ytick.gridOn': grid,
          'ytick.major.width':5,
          'xtick.major.width':5,
          'ytick.major.size': 5,
          'xtick.major.size': 5,
          'mathtext.fontset': 'stixsans',
          'mathtext.sf': 'sans',
          'legend.frameon': True,
          'legend.facecolor': '#FFEDCE',
          'figure.dpi': 150,
           'xtick.color': 'k',
           'ytick.color': 'k'}
    plt.rc('text.latex', preamble=r'\usepackage{sfmath}')
    plt.rc('mathtext', fontset='stixsans', sf='sans')
    sns.set_style('darkgrid', rc=rc)

