using PlotlyJS

# Set default values for all plots
nlaStyle = let
    axis = attr(showgrid=true, gridcolor="#E5E5E5",
        linewidth=1.0,
        titlefont_color="#555555", titlefont_size=18,
        linecolor="black", mirror=true, zeroline=false,
        ticks="inside")
    layout = Layout(font_size=16, xaxis=axis, yaxis=axis,
        titlefont_size=18, width=500, height=300)

    colors = Cycler([
            "#E24A33", "#348ABD", "#988ED5", "#777777",
            "#FBC15E", "#8EBA42", "#FFB5B8"
    ])
    gta = attr(
        marker_line_width=0.5, marker_line_color="#348ABD",
        marker_color=colors, marker_size=10
    )
    Style(layout=layout, global_trace=gta)
end

use_style!(nlaStyle)

function plotToPDF(p,fileName)
    savefig(p,"$fileName.html")
    cmd = `/usr/libexec/cups/filter/xhtmltopdf 0 0 0 0 0 $fileName.html`
    run(pipeline(cmd, stdout="$fileName.pdf"))
    cmd = `pdfcrop $fileName.pdf`
    run(cmd)
end
