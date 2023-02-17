import numpy as np
import pylab as plt
import pandas as pd
import sys

if not ("-f" in sys.argv and "-o" in sys.argv and len(sys.argv) == 5):
    print("USAGE:\n"+sys.argv[0]+" -f 2d.xvg -o out_path")
    exit()
    
for i in range(1, len(sys.argv)) :
    if sys.argv[i] == "-f":
        i_file = sys.argv[i+1]
    if sys.argv[i] == "-o":
        o_file = sys.argv[i+1]

data = np.loadtxt(i_file, skiprows=17)

df = {'e1':data[:, 0], 'e2':data[:, 1]}
df = pd.DataFrame(data=df)

import plotly.express as px

fig = px.density_contour(df, y="e1", x="e2", width=1000, height=800)
fig.update_yaxes(range=[-15, 20], title="projection on eigenvector 1 (nm^2)")
fig.update_xaxes(range=[-10, 15], title="projection on eigenvector 2 (nm^2)")
fig.update_traces(contours_coloring="fill", contours_showlabels = True, colorscale='Hot_r')
fig.update_layout(plot_bgcolor='white')

fig.write_html(o_file+'.html')
