import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import scipy.optimize as optimize
from pyparsing import results
from scipy.optimize import least_squares
import math
from dataclasses import dataclass, field

@dataclass
class car:
    total_mass: float #in Kgs
    sprung_mass: float #in Kgs
    CG_height: float
    track_width: float
    ride_height: float
    upper_inner_R: tuple[float, float]
    lower_inner_R: tuple[float, float]
    length_upper: float
    length_lower: float
    length_upright: float
    spring_stiffness: float
    sway_bar_stiffness: float
    upper_inner_L: tuple[float, float] = field(init=False)
    lower_inner_L: tuple[float, float] = field(init=False)

    #all units otherwise mentioned are in mm
    def __post_init__(self):
        self.upper_inner_L = (-self.upper_inner_R[0], self.upper_inner_R[1])
        self.lower_inner_L = (-self.lower_inner_R[0], self.lower_inner_R[1])

def suspension_eqns(vars, upper_inner, lower_inner, L_upper, L_lower, L_upright, WHEELCENTER):
    x_Uo, y_Uo, x_Lo, y_Lo = vars

    eq1 = ((x_Uo - upper_inner[0])**2 + (y_Uo - upper_inner[1])**2)**0.5 - L_upper

    eq2 = ((x_Lo - lower_inner[0])**2 + (y_Lo - lower_inner[1])**2)**0.5 - L_lower

    eq3 = ((x_Uo - x_Lo)**2 + (y_Uo - y_Lo)**2)**0.5 - L_upright

    eq4 = (y_Uo + y_Lo)/2 - WHEELCENTER[1]

    # print(eq1, eq2, eq3, eq4)
    return [eq1, eq2, eq3, eq4]

def solve_position(guess, u_in, l_in, L_u, L_l, L_up, WHEELCENTER, solver):
    if solver == 0:
        solution, info, ier, msg = optimize.fsolve(suspension_eqns, x0=guess, args=(u_in, l_in, L_u, L_l, L_up, WHEELCENTER), full_output=True)
        # print(ier)
        return solution
    elif solver == 1:
        solution = least_squares(suspension_eqns, x0=guess, args=(u_in, l_in, L_u, L_l, L_up, WHEELCENTER))
        return solution.x

def get_camber(x_Uo, y_Uo, x_Lo, y_Lo):
    upright_vec = np.array([x_Uo - x_Lo, y_Uo - y_Lo])
    camber = np.sign(upright_vec[0])*math.degrees(math.acos(upright_vec[1]/np.linalg.norm(upright_vec)))
    return camber

def intersection(p1, p2, p3, p4):
    x1, y1= p1
    x2, y2 = p2
    x3, y3 = p3
    x4, y4 = p4

    denominator = (x1 - x2)*(y3 - y4) - (y1 - y2)*(x3 - x4)

    if abs(denominator) < 1e-6:
        return None

    px = ((x1*y2 - y1*x2)*(x3-x4) - (x1 - x2)*(x3*y4 - y3*x4))/denominator
    py = ((x1*y2 - y1*x2)*(y3-y4) - (y1 - y2)*(x3*y4 - y3*x4))/denominator

    return px, py

def y_at_x(p1, p2, x):
    x1, y1 = p1
    x2, y2 = p2

    if abs(x2 - x1) < 1e-6:
        return None

    slope = (y2 - y1)/(x2 - x1)
    return y1 + slope * (x - x1)

def simulate(car, span):
    guess_R = None
    guess_L = None
    dy_vals_R = []
    dy_vals_L = []
    camber_vals_R = []
    camber_vals_L = []
    rc_vals_R = []
    rc_vals_L = []
    for dy in range(-span, span+1):

        dy_R = dy
        dy_L = -dy

        wc_R = (car.track_width / 2, car.ride_height + dy_R)
        wc_L = (-car.track_width / 2, car.ride_height + dy_L)

        if guess_R is None:
            guess_R = [
                # SUBSTITUTE RIGHT SIDE SUSPENSION GUESSES
                599.74259, # upper outer x
                399.99978, # upper outer y
                600.25741, # lower outer x
                100.00022  # lower outer y
            ]
        if guess_L is None:
            guess_L = [
                # SUBSTITUTE LEFT SIDE SUSPENSION GUESSES
                -599.74259, # upper outer x
                399.99978,  # upper outer y
                -600.25741, # lower outer x
                100.00022   # lower outer y
            ]

        sol_R = solve_position(guess_R, car.upper_inner_R, car.lower_inner_R, car.length_upper, car.length_lower, car.length_upright, wc_R, 0)
        upper_outer_R = (sol_R[0], sol_R[1])
        lower_outer_R = (sol_R[2], sol_R[3])

        sol_L = solve_position(guess_L, car.upper_inner_L, car.lower_inner_L, car.length_upper, car.length_lower, car.length_upright, wc_L, 0)
        upper_outer_L = (sol_L[0], sol_L[1])
        lower_outer_L = (sol_L[2], sol_L[3])

        #Camber Simulation
        camber_R = get_camber(upper_outer_R[0], upper_outer_R[1], lower_outer_R[0], lower_outer_R[1])
        dy_vals_R.append(dy_R)
        camber_vals_R.append(camber_R)

        camber_L = get_camber(upper_outer_L[0], upper_outer_L[1], lower_outer_L[0], lower_outer_L[1])
        dy_vals_L.append(dy_L)
        camber_vals_L.append(camber_L)

        #Roll Center Simulation
        contactPatch_R = (wc_R[0], 0)
        contactPatch_L = (wc_L[0], 0)
        IC_R = intersection(car.upper_inner_R, (sol_R[0], sol_R[1]), car.lower_inner_R, (sol_R[2], sol_R[3]))
        if IC_R:
            rc_height_R = y_at_x(IC_R, contactPatch_R, 0)
            rc_vals_R.append(rc_height_R)
        IC_L = intersection(car.upper_inner_L, (sol_L[0], sol_L[1]), car.lower_inner_L, (sol_L[2], sol_L[3]))
        if IC_L:
            rc_height_L = y_at_x(IC_L, contactPatch_L, 0)
            rc_vals_L.append(rc_height_L)

        guess_R = sol_R
        guess_L = sol_L

        if dy ==0:
            static_L = sol_L
            static_R = sol_R
            static_rc_L = rc_height_L
            static_rc_R = rc_height_R
            print(rc_height_R, rc_height_L)
            print(IC_R, IC_L)
            print(car.upper_inner_R, (sol_R[0], sol_R[1]), car.lower_inner_R, (sol_R[2], sol_R[3]))

    return (dy_vals_R, camber_vals_R, rc_vals_R, sol_R), (dy_vals_L, camber_vals_L, rc_vals_L, sol_L), (static_L, static_R, static_rc_L, static_rc_R)

def main():
    mycar = car(
        #EDIT THESE VALUES TO YOUR USE CASE
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
    )
    results_R, results_L, statics= simulate(mycar, 20)
    dy_vals_R, camber_vals_R, rc_vals_R, sol_R = results_R
    dy_vals_L, camber_vals_L, rc_vals_L, sol_L = results_L
    static_L, static_R, static_rc_L, static_rc_R = statics

    # Camber Curve
    plt.subplot(3, 2, 1)
    plt.plot(dy_vals_L, camber_vals_L, color="b")
    plt.xlabel("Suspension travel (mm)")
    plt.ylabel("Camber (deg)")
    plt.grid(True)

    plt.subplot(3, 2, 2)
    plt.plot(dy_vals_R, camber_vals_R, color="r")
    plt.xlabel("Suspension travel (mm)")
    plt.ylabel("Camber (deg)")
    plt.grid(True)

    # Roll Center Curve
    plt.subplot(3, 2, 3)
    plt.plot(dy_vals_L, rc_vals_L, color="b")
    plt.xlabel("Suspension travel (mm)")
    plt.ylabel("Roll center height (mm)")
    plt.grid(True)

    plt.subplot(3, 2, 4)
    plt.plot(dy_vals_R, rc_vals_R, color="r")
    plt.xlabel("Suspension travel (mm)")
    plt.ylabel("Roll center height (mm)")
    plt.grid(True)

    #Visual representation of final suspension state
    plt.subplot(3, 2, 5)
    # Left side
    plt.scatter(mycar.upper_inner_L[0], mycar.upper_inner_L[1], label="upper inner")
    plt.scatter(mycar.lower_inner_L[0], mycar.lower_inner_L[1], label="lower inner")
    plt.scatter(static_L[0], static_L[1], label="upper outer")
    plt.scatter(static_L[2], static_L[3], label="lower outer")
    plt.scatter((static_L[0] + static_L[2]) / 2, (static_L[1] + static_L[3]) / 2, label="wheel center")
    plt.scatter(-mycar.track_width/2, static_rc_L, label="Roll Center")
    plt.plot((mycar.upper_inner_L[0], static_L[0]), (mycar.upper_inner_L[1], static_L[1]), label="upper arm")
    plt.plot((mycar.lower_inner_L[0], static_L[2]), (mycar.lower_inner_L[1], static_L[3]), label="lower arm")
    plt.plot((static_L[0], static_L[2]), (static_L[1], static_L[3]), label="upright")
    # Right side
    plt.scatter(mycar.upper_inner_R[0], mycar.upper_inner_R[1], label="upper inner")
    plt.scatter(mycar.lower_inner_R[0], mycar.lower_inner_R[1], label="lower inner")
    plt.scatter(static_R[0], static_R[1], label="upper outer")
    plt.scatter(static_R[2], static_R[3], label="lower outer")
    plt.scatter((static_R[0] + static_R[2]) / 2, (static_R[1] + static_R[3]) / 2, label="wheel center")
    plt.scatter(mycar.track_width / 2, static_rc_R, label="Roll Center")
    plt.plot((mycar.upper_inner_R[0], static_R[0]), (mycar.upper_inner_R[1], static_R[1]), label="upper arm")
    plt.plot((mycar.lower_inner_R[0], static_R[2]), (mycar.lower_inner_R[1], static_R[3]), label="lower arm")
    plt.plot((static_R[0], static_R[2]), (static_R[1], static_R[3]), label="upright")
    plt.grid(True)

    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.98, top=0.98, wspace=0.2, hspace=0.36)
    # plt.tight_layout()
    plt.show()

main()
