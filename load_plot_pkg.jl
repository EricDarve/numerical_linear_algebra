using PlotlyJS

function set_plotly_template()
    axis = attr(showgrid=true, gridcolor="#E5E5E5",
        linewidth=1.0,
        title_font_color="#555555", title_font_size=14,
        linecolor="black", showline=true, mirror=true, zeroline=false,
        ticks="inside")
    layout=Layout(font_size=16, xaxis=axis, yaxis=axis,
        titlefont_size=18, width=500, height=300)

    layout[:colorway]=[
            "#E24A33", "#348ABD", "#988ED5", "#777777",
            "#FBC15E", "#8EBA42", "#FFB5B8"]
    
    t = Template(layout=layout)
    t.data[:scatter]  = [attr(marker_line_width=0.5,
                              marker_line_color="#348ABD",
                              marker_size=10)]
    return t
end

fig_template = set_plotly_template()

# CHANGE: add the template under the "personal" key
templates.personal = fig_template

# CHANGE: set the default template to the "personal" template -- must match key from previous step
templates.default = "personal"

function plotToPDF(p,fileName)
    savefig(p,"$fileName.html")
    cmd = `/usr/libexec/cups/filter/xhtmltopdf 0 0 0 0 0 $fileName.html`
    run(pipeline(cmd, stdout="$fileName.pdf"))
    cmd = `pdfcrop $fileName.pdf`
    run(cmd)
end