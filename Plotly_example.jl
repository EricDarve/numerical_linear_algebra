import Plotly
line1  = Plotly.scatter(x=1:n, y=dA[:,1], mode="lines", name="Pivoted LU")
line2  = Plotly.scatter(x=1:n, y=dA[:,2], mode="lines", name="Rook LU")
layout = Plotly.Layout(xaxis_title="Diagonal entry", yaxis_title="Magnitude", yaxis_type = "log", 
                                                                             margin=Dict("l"=>50))
#println(JSON.json(layout, 2))
#margin=Dict("pad" => 0)
#yaxis_exponentformat
#margin=Dict() 
Plotly.plot([line1,line2], layout)
