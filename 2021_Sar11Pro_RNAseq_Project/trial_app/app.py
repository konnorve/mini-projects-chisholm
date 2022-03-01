import dash
from dash import dcc, html, callback_context

app = dash.Dash(__name__)
server = app.server

dropdown_options = [{'label':'', 'value':'none'},
                            {'label':'Control', 'value':'control'},
                            {'label':'Treatment', 'value':'treatment'},]

app.layout = html.Div([
    html.H4('Hello World', style={'text-align': 'center'}),
    dcc.Dropdown(options=dropdown_options, value='none')
], style={'display': 'grid', 'grid-template-columns': '1fr 2fr'})

if __name__ == "__main__":
    app.run_server(debug=True, host='127.0.0.1', port='8000')