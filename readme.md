# Double Wishbone Suspension Simulator
A python based simulator for symmetric double wishbone suspensions to simulate camber curve and roll center height through a range of suspension travel.

### What is this?
This program helps users determine a suitable double wishbone suspension geometry for their use case by simulating the state of the suspension at static ride height and plotting a camber curve and a roll center height curve through a specified suspension travel range. 

### How do I run it?
This program uses various python modules that must first be installed to your python interpreter before execution. The modules used are scipy, matplotlib and numpy. 

##### To install, the dependencies, run the following command in the same directory that contains main.py and requirements.txt:

    pip install -r requirements.txt

Once the dependencies are installed, download the main.py file and run it. This should run the simulator with the default example suspension geometry values.

### How do I use it?
##### *Note: This simulator uses a 2D Cartesian system for coordinates of hardpoints with x = 0 as the car centerline and y = 0 as the road surface. Only right side hardpoints are required as the left side is computed automatically by mirroring.*

This simulator takes a number of inputs to provide an output. They are as shown below:

        total_mass: float #Total mass of the vehicle in kgs
        sprung_mass: float #Sprung mass of the vehicle in kgs
        CG_height: float #Height of center of gravity of the vehicle
    
        track_width: float #Track width of vehicle (distance between contact patches on opposite tyres)
        ride_height: float #Static ride height of the vehicle
    
        upper_inner_R: tuple[float, float] #position of the upper inner pivot of the suspension
        lower_inner_R: tuple[float, float] #position of the lower inner pivot of the suspension
    
        length_upper: float #length of the upper control arm
        length_lower: float	#length of lower control arm
        length_upright: float #length of upright
    
        spring_stiffness: float #stiffness of the suspension spring in N/mm
        sway_bar_stiffness: float #stiffness of the sway bar in N/mm

*Note: all fields above are in mm unless mentioned otherwise.*

The default example values are as follows:

	    total_mass=300,
        sprung_mass=240,
        CG_height=300,

        track_width=1200,
        ride_height=250,

        upper_inner_R=(260, 340),
        lower_inner_R=(260, 180),

        length_upper=345,
        length_lower=350,
        length_upright=300,

        spring_stiffness=35,
        sway_bar_stiffness=15

To change the simulation parameters, change the values for each field under the comment `#EDIT THESE VALUES TO YOUR USE CASE`

Once the simulation parameters are changed, an appropriate guess for the suspension's outer pivot positions must be provided. To arrive at an accurate guess, plot the hardpoints into a graphical calculator like desmos and find the approximate positions for the outer pivots. Then, substitute the guesses under the comment `#SUBSTITUTE <SIDE> SIDE SUSPENSION GUESSES`

### What does it output?
The simulator uses matplotlib to neatly plot and display the camber curves and roll center height curves for the left and right side suspensions individually in subplots. It also displays a visualization of the suspension geometry at static ride height along with markers for the roll center height at static ride height.

![Example Simulator Output](example%20output.png)

### Limitations and assumptions:
- Assumes suspension is symmetric
- Simulation is in 2D
- Assumes steady state
- No tire model is considered
- Static center of gravity is assumed to be at the vehicle center line
- Requires accurate initial guesses for fsolve to converge correctly
- No easy to use UI to set suspension geometry to be tested


