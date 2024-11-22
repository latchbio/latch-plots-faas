import plotly.graph_objects as go


def graphpad_inspired_theme() -> go.layout.Template:
    border_width = 3
    border_color = "#000000"
    font_color = "#000000"

    axis_common = {
        "title": {"font": {"size": 16, "color": font_color}},
        "showline": True,
        "linecolor": border_color,
        "linewidth": border_width,
        "ticks": "outside",
        "ticklen": 5,
        "tickwidth": border_width,
        "tickcolor": border_color,
    }

    return go.layout.Template(
        layout={
            "title": {
                "font": {
                    "family": "InterVariable, Inter, sans-serif",
                    "size": 24,
                    "weight": 700,
                    "color": font_color,
                }
            },
            "font": {
                "family": "InterVariable, Inter, sans-serif",
                "weight": 700,
                "size": 16,
                "color": font_color,
            },
            "colorway": ["#0634FF", "#FF2700", "#04C703", "#BD38E9", "#FF9301"],
            "xaxis": axis_common,
            "yaxis": axis_common,
            "bargap": 0.30,
            "bargroupgap": 0.30,
        },
        data={
            "bar": [
                go.Bar(marker={"line": {"width": border_width, "color": border_color}})
            ]
        },
    )
