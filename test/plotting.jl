using Plots



"""
Some default parameters for plotting
"""

FONT = "Computer Modern" 
TITLESIZE = 10
TITLEFONT = Plots.font(TITLESIZE, FONT)
GUIDESIZE = 10
GUIDEFONT = Plots.font(GUIDESIZE, FONT)
TICKSIZE = 10
TICKFONT = Plots.font(TICKSIZE, FONT)
LEGENDSIZE = 10
LEGENDFONT = Plots.font(LEGENDSIZE, FONT)

GRID_DEFAULT = false
SIZE_DEFAULT = (350, 350)


function my_plot(title = ""; 
    titlefont = TITLEFONT,
    guidfont = GUIDEFONT, 
    tickfont = TICKFONT, 
    legendfont = LEGENDFONT,
    grid = GRID_DEFAULT,
    size = SIZE_DEFAULT)
    return plot(title = title, 
    titlefont = TITLEFONT,
    guidfont = GUIDEFONT,
    tickfont = TICKFONT,
    legendfont = LEGENDFONT, 
    grid = GRID_DEFAULT,
    size = SIZE_DEFAULT)
end 



