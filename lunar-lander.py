import matplotlib.pyplot as plt
import matplotlib.widgets as widgets
import matplotlib.patches as patches
from matplotlib.path import Path

# Functions used by the program
def closeCallback(event):
    plt.close('all') # Close all open figure windows
    exit()

def thrustCallback(val):
    global thrust
    thrust = val

def drawLander(x, y):
    verts = [
        (x-3.,y),
        (x-1.5,y+2),
        (x-3.,y+2),
        (x-3.,y+10.),
        (x,y+14.),
        (x+3.,y+10.),
        (x+3.,y+2),
        (x+1.5,y+2),
        (x+3.,y),
        (x-3.,y)
    ]
    codes = [ Path.MOVETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.LINETO, Path.CLOSEPOLY ]
    path = Path(verts, codes)
    patch = patches.PathPatch(path, facecolor='gray')
    ax.add_patch(patch)

# Setting system parameters for the game
dt = 1/30 # timestep for simulation
h0 = 100 # starting height for lander in meters
g = -1.0 # setting the gravitational acceleration in m/s^2
maxthrust = 3.0 # setting the maximum engine acceleration in m/s^2
crashspeed = 3.0 # setting the speed at which the landing is considered a crash in m/s

# Gui axes
fig = plt.figure(figsize=(7, 7))
ax = plt.axes([0.1, 0.2, 0.6, 0.65]) # Axes of main plot
close_ax = plt.axes([0.8, 0.775, 0.1, 0.1]) # Axes of close button box
thrust_ax = plt.axes([0.15, 0.075, 0.7, 0.05]) # Axes of thrust slider
fuel_ax = plt.axes([0.775, 0.2, 0.15, 0.5]) # Axes of the fuel indicator
# gui widgets
closeHandle = widgets.Button(close_ax, 'Close') # Close button
closeHandle.on_clicked(closeCallback)
thrustHandle = widgets.Slider(thrust_ax, 'Thrust', 0, maxthrust, valinit=0)
thrustHandle.on_changed(thrustCallback)

# Initialising variables
t = 0 # setting the starting time of the simulstion to t=0 seconds
x = h0 # setting initial position of the lander
v = 0 # setting the initial speed of the lander in m/s
fuel = 100 # starting fuel
thrust = 0

# Main logic loop
while x > 0:
    if fuel > 0:
        a = g + thrust
    else:
        a = g
    dx = v*dt
    dv = a*dt
    x += dx
    v += dv
    t += dt
    fuel -= thrust/4
    ax.cla()
    ax.set_ylim(0,h0+h0/5)
    ax.set_xlim(-h0/2,h0/2)
    if -v >= crashspeed:
        ax.set_title('Velocity: {0:.2f} m/s   Altitude: {1:.2f} m'.format(v,x), fontsize = 15, color='red')
    else:
        ax.set_title('Velocity: {0:.2f} m/s   Altitude: {1:.2f} m'.format(v,x), fontsize = 15)
    drawLander(0,x)
    fuel_ax.cla()
    fuel_ax.set_title("Fuel")
    fuel_ax.set_ylim(0,100)
    fuel_ax.set_xlim(0,1)
    fuel_ax.bar([0.5], [fuel], 1, color='orange')
    plt.pause(dt)

# Was it a crash landing?
if -v > crashspeed:
    ax.text(-18, 30, "You crashed!", fontsize=18)
    plt.pause(dt)
else:
    ax.text(-32, 30, "Another happy landing", fontsize=18)
    plt.pause(dt)

# Show the graph
plt.show()
