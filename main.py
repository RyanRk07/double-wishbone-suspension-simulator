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
    total_sprung_mass: float #in Kgs
    front_sprung_mass: float = field(init=False) #in Kgs
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
    mount_fraction: float
    shock_mount_upper_R: tuple[float, float]
    shock_mount_upper_L: tuple[float, float] = field(init=False)
    upper_inner_L: tuple[float, float] = field(init=False)
    lower_inner_L: tuple[float, float] = field(init=False)

    #all units otherwise mentioned are in mm
    def __post_init__(self):
        self.upper_inner_L = (-self.upper_inner_R[0], self.upper_inner_R[1])
        self.lower_inner_L = (-self.lower_inner_R[0], self.lower_inner_R[1])
        self.shock_mount_upper_L = (-self.shock_mount_upper_R[0], self.shock_mount_upper_R[1])
        self.front_sprung_mass = self.total_sprung_mass/2

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

def distance(p1, p2):
    x1, y1 = p1
    x2, y2 = p2
    return math.sqrt((x1 - x2)**2 + (y1 - y2)**2)

def shock_mount_pos(lower_inner, lower_outer, mount_fraction):
    y = lower_outer[1]*(mount_fraction) + (lower_inner[1]*(1 - mount_fraction))
    x = lower_outer[0]*(mount_fraction) + (lower_inner[0]*(1 - mount_fraction))
    return (x, y)

def calculate_installation_angle(lower_inner, lower_outer, shock_lower, shock_upper):
    arm = np.array(lower_outer) - np.array(lower_inner)
    shock = np.array(shock_upper) - np.array(shock_lower)
    norm_arm = np.linalg.norm(arm)
    norm_shock = np.linalg.norm(shock)
    dotprod = np.dot(arm, shock)
    cosine = np.clip(dotprod / (norm_arm*norm_shock), -1.0, 1.0)
    angle = np.degrees(np.arccos(cosine))
    if angle > 90:
        angle = 180 - angle
    return angle

def simulate(car, a):
    guess_L = None
    guess_R = None
    dy_vals_R = []
    dy_vals_L = []
    camber_vals_R = []
    camber_vals_L = []
    rc_vals_R = []
    rc_vals_L = []
    instant_radii_L = []
    instant_radii_R = []
    shock_mount_pos_R = []
    shock_mount_pos_L = []
    spring_lengths_L = []
    spring_lengths_R = []
    installation_angles_R = []
    installation_angles_L = []
    wheel_rates_R = []
    wheel_rates_L = []
    instant_loads_R = []
    instant_loads_L = []

    cornering_acceleration = a * 9.18
    lateral_force = car.front_sprung_mass * cornering_acceleration
    load_transferred = (lateral_force * car.CG_height)/car.track_width
    print(load_transferred)

    target_load = load_transferred
    current_load_R = 0
    current_load_L = 0
    step = 0.1 * np.sign(a)
    dy_R = 0
    dy_L = 0

    while current_load_R < target_load:

        wc_R = (car.track_width / 2, car.ride_height + dy_R)
        contactPatch_R = (wc_R[0], 0)

        if guess_R is None:
            guess_R = [
                # SUBSTITUTE RIGHT SIDE SUSPENSION GUESSES
                599.74259, # upper outer x
                399.99978, # upper outer y
                600.25741, # lower outer x
                100.00022  # lower outer y
            ]

        sol_R = solve_position(guess_R, car.upper_inner_R, car.lower_inner_R, car.length_upper, car.length_lower, car.length_upright, wc_R, 0)
        upper_outer_R = (sol_R[0], sol_R[1])
        lower_outer_R = (sol_R[2], sol_R[3])
        guess_R = sol_R

        camber_R = get_camber(upper_outer_R[0], upper_outer_R[1], lower_outer_R[0], lower_outer_R[1])
        camber_vals_R.append(camber_R)

        IC_R = intersection(car.upper_inner_R, (sol_R[0], sol_R[1]), car.lower_inner_R, (sol_R[2], sol_R[3]))
        if IC_R:
            rc_height_R = y_at_x(IC_R, contactPatch_R, 0)
            rc_vals_R.append(rc_height_R)
            instant_radii_R.append(distance(IC_R, wc_R))

        shock_mount_lower_R = shock_mount_pos(car.lower_inner_R, lower_outer_R, car.mount_fraction)
        shock_mount_pos_R.append(shock_mount_lower_R)

        spring_lengths_R.append(distance(car.shock_mount_upper_R, shock_mount_lower_R))

        installation_angle = calculate_installation_angle(car.lower_inner_R, lower_outer_R, shock_mount_lower_R, car.shock_mount_upper_R)
        installation_angles_R.append(installation_angle)

        if dy_R == 0:
            static_R = sol_R
            static_rc_R = rc_height_R
            static_shock_mount_R = shock_mount_lower_R
            static_spring_R = distance(static_shock_mount_R, car.shock_mount_upper_R)

        if dy_R == 0:
            motion_ratio = (distance(car.shock_mount_upper_R, static_shock_mount_R)/distance(car.shock_mount_upper_R, lower_outer_R))*np.sin(installation_angle)
        else:
            motion_ratio = (distance(car.shock_mount_upper_R, shock_mount_lower_R) - static_spring_R)/dy_R

        wheel_rate = car.spring_stiffness * (motion_ratio**2)
        wheel_rates_R.append(wheel_rate)
        instant_loads_R.append(current_load_R)
        current_load_R += wheel_rate*step
        dy_vals_R.append(dy_R)
        dy_R += step


    while current_load_L < target_load:

        wc_L = (-car.track_width / 2, car.ride_height + dy_L)
        contactPatch_L = (wc_L[0], 0)

        if guess_L is None:
            guess_L = [
                # SUBSTITUTE LEFT SIDE SUSPENSION GUESSES
                -599.74259, # upper outer x
                399.99978,  # upper outer y
                -600.25741, # lower outer x
                100.00022   # lower outer y
            ]

        sol_L = solve_position(guess_L, car.upper_inner_L, car.lower_inner_L, car.length_upper, car.length_lower, car.length_upright, wc_L, 0)
        upper_outer_L = (sol_L[0], sol_L[1])
        lower_outer_L = (sol_L[2], sol_L[3])
        guess_L = sol_L

        camber_L = get_camber(upper_outer_L[0], upper_outer_L[1], lower_outer_L[0], lower_outer_L[1])
        camber_vals_L.append(camber_L)

        IC_L = intersection(car.upper_inner_L, (sol_L[0], sol_L[1]), car.lower_inner_L, (sol_L[2], sol_L[3]))
        if IC_L:
            rc_height_L = y_at_x(IC_L, contactPatch_L, 0)
            rc_vals_L.append(rc_height_L)
            instant_radii_L.append(distance(IC_L, wc_L))

        shock_mount_lower_L = shock_mount_pos(car.lower_inner_L, lower_outer_L, car.mount_fraction)
        shock_mount_pos_L.append(shock_mount_lower_L)

        spring_lengths_L.append(distance(car.shock_mount_upper_L, shock_mount_lower_L))

        installation_angle = calculate_installation_angle(car.lower_inner_L, lower_outer_L, shock_mount_lower_L,
                                                          car.shock_mount_upper_L)
        installation_angles_L.append(installation_angle)

        if dy_L == 0:
            static_L = sol_L
            static_rc_L = rc_height_L
            static_shock_mount_L = shock_mount_lower_L
            static_spring_L = distance(static_shock_mount_L, car.shock_mount_upper_L)

        if dy_L == 0:
            motion_ratio = (distance(car.shock_mount_upper_L, static_shock_mount_L) / distance(car.shock_mount_upper_L, lower_outer_L)) * np.sin(installation_angle)
        else:
            motion_ratio = (distance(car.shock_mount_upper_L, shock_mount_lower_L) - static_spring_L) / dy_L

        wheel_rate = car.spring_stiffness * (motion_ratio ** 2)
        wheel_rates_L.append(wheel_rate)
        instant_loads_L.append(current_load_L)
        current_load_L += wheel_rate * step
        dy_vals_L.append(dy_L)
        dy_L -= step

    gforces_R = (np.array(instant_loads_R)*car.track_width)/(9.18*car.CG_height*car.front_sprung_mass)
    gforces_L = (np.array(instant_loads_L)*car.track_width)/(9.18*car.CG_height*car.front_sprung_mass)
    limiter = min(len(gforces_L), len(gforces_R))
    roll_angles = np.atan((abs(np.array(dy_vals_L[:limiter])) + abs(np.array(dy_vals_R[:limiter]))) / car.track_width)
    return (dy_vals_R, camber_vals_R, rc_vals_R, sol_R), (dy_vals_L, camber_vals_L, rc_vals_L, sol_L), (static_L, static_R, static_rc_L, static_rc_R, static_spring_L, static_spring_R), (instant_radii_L, instant_radii_R), (spring_lengths_L, spring_lengths_R), (installation_angles_L, installation_angles_R), (wheel_rates_L, wheel_rates_R), (instant_loads_L, instant_loads_R), roll_angles, (gforces_L, gforces_R)

def main():
    mycar = car(
        #EDIT THESE VALUES TO YOUR USE CASE
        total_sprung_mass=200,
        CG_height=300,

        track_width=1200,
        ride_height=250,

        upper_inner_R=(260, 340),
        lower_inner_R=(260, 180),

        length_upper=345,
        length_lower=350,
        length_upright=300,

        spring_stiffness=35,
        sway_bar_stiffness=15,
        mount_fraction=0.6,
        shock_mount_upper_R=(464, 282)
    )
    results_R, results_L, statics, instant_radii, spring_lengths, installation_angles, wheel_rates, instant_loads, roll_angles, gforces = simulate(mycar, 1.0)
    gforces_L, gforces_R = gforces
    instant_loads_L, instant_loads_R = instant_loads
    installation_angles_L, installation_angles_R = installation_angles
    instant_radii_L, instant_radii_R = instant_radii
    spring_lengths_L, spring_lengths_R = spring_lengths
    dy_vals_R, camber_vals_R, rc_vals_R, sol_R = results_R
    dy_vals_L, camber_vals_L, rc_vals_L, sol_L = results_L
    static_L, static_R, static_rc_L, static_rc_R, static_spring_L, static_spring_R = statics
    wheel_rates_L, wheel_rates_R = wheel_rates

    print(len(dy_vals_R), len(dy_vals_L), dy_vals_R[-1], dy_vals_L[-1], dy_vals_R[0], dy_vals_L[0])
    print(len(wheel_rates_R), len(wheel_rates_L), wheel_rates_R[-1], wheel_rates_L[-1], wheel_rates_R[0], wheel_rates_L[0])
    print(instant_loads_L[-1], instant_loads_R[-1])


    # Camber Curve
    plt.subplot(4, 2, 1)
    plt.plot(gforces_L, camber_vals_L, color="b", alpha = 0.5)
    plt.plot(gforces_R, camber_vals_R, color="r", alpha = 0.5)
    plt.text(10, 0, f"camber gain = {camber_vals_R[-1]/dy_vals_R[-1]:.2f}")
    plt.xlabel("g forces (g)")
    plt.ylabel("Camber (deg)")
    plt.grid(True)


    # Roll Center Curve
    plt.subplot(4, 2, 2)
    plt.plot(gforces_L, rc_vals_L, color="b", alpha = 0.5)
    plt.plot(gforces_R, rc_vals_R, color="r", alpha = 0.5)
    plt.xlabel("g forces (g)")
    plt.ylabel("Roll center height (mm)")
    plt.grid(True)


    #Visual representation of final suspension state
    plt.subplot(4, 2, 3)
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

    plt.subplot(4, 2, 4)
    plt.plot(gforces_R, instant_radii_R, color="r", alpha = 0.5)
    plt.plot(gforces_L, instant_radii_L, color="b", alpha = 0.5)
    plt.xlabel("g forces (g)")
    plt.ylabel("Instant radius (mm)")
    plt.grid(True)

    plt.subplot(4, 2, 5)
    plt.plot(gforces_L, wheel_rates_L, color="b", alpha = 0.5)
    plt.plot(gforces_R, wheel_rates_R, color="r", alpha = 0.5)
    plt.xlabel("g forces (g)")
    plt.ylabel("Wheel rate (N/mm)")
    plt.grid(True)

    plt.subplot(4, 2, 6)
    plt.plot(gforces_R, installation_angles_R, color="r", alpha = 0.5)
    plt.plot(gforces_L, installation_angles_L, color="b", alpha = 0.5)
    plt.xlabel("g forces (g)")
    plt.ylabel("Installation angle (deg)")
    plt.grid(True)

    plt.subplot(4, 2, 7)
    plt.plot(gforces_L, instant_loads_L, color="b", alpha = 0.5)
    plt.plot(gforces_R, instant_loads_R, color="r", alpha = 0.5)
    plt.xlabel("g forces (g)")
    plt.ylabel("Change in load (N)")
    plt.grid(True)

    plt.subplot(4, 2, 8)
    if len(gforces_R) < len(gforces_L):
        plt.plot(gforces_R, roll_angles, color="g", alpha = 0.5)
    else:
        plt.plot(gforces_L, roll_angles, color="g", alpha = 0.5)
    plt.xlabel("g forces (g)")
    plt.ylabel("Body roll angle (deg)")
    plt.grid(True)

    plt.subplots_adjust(left=0.06, bottom=0.06, right=0.98, top=0.98, wspace=0.2, hspace=0.36)
    # plt.tight_layout()
    plt.show()

main()
