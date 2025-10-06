from vpython import *
import math
import datetime

# --- Scene Setup ---
scene.title = "Epoch-Aligned Elliptical Solar System"
scene.range = 35  # Set the initial viewing range to see Neptune/Pluto (AU)
scene.autoscale = False
# Set scene to view the Ecliptic plane (orbits) in the X-Y plane (horizontal)
scene.forward = vector(0, -1, 0) # Camera looks along the negative Y-axis
scene.up = vector(0, 0, 1)      # Z-axis is 'up' (North Ecliptic Pole)
scene.background = color.black

# --- Constants ---
G_DAYS = 365.25 # Days in a year for period calculation
J2000_JD = 2451545.0

# --- Time Control ---
sim_days = 0.0
time_scale = 100000.0 # Default speed: 100,000x real time
last_time = datetime.datetime.now()

# --- Utility Functions ---

def jd_of_date(dt):
    """Calculate Julian Day for a Python datetime object."""
    # Days since 1970-01-01 00:00:00 UTC (The Unix Epoch)
    seconds_since_epoch = dt.timestamp()
    # JD of Unix Epoch is 2440587.5
    return seconds_since_epoch / 86400.0 + 2440587.5

def days_since_j2000(dt):
    """Calculate days since J2000 epoch (JD - J2000_JD)."""
    return jd_of_date(dt) - J2000_JD

def deg2rad(d):
    """Convert degrees to radians."""
    return d * math.pi / 180

def norm_deg(d):
    """Normalize angle to [0, 360)."""
    d %= 360
    return d if d >= 0 else d + 360

# --- Kepler Solver (Newton-Raphson) ---

def solve_kepler(M, e):
    """Solve Kepler's equation M = E - e*sin(E) for Eccentric Anomaly E."""
    M = math.atan2(math.sin(M), math.cos(M))
    E = M + e * math.sin(M) * (1.0 + e * math.cos(M))
    
    for _ in range(12):
        f = E - e * math.sin(E) - M
        fp = 1 - e * math.cos(E)
        dE = f / fp
        E -= dE
        if abs(dE) < 1e-12:
            break
    return E

# --- Orbital Mechanics ---

def propagate_m_deg(M0_deg, epochJD, targetJD, a_AU):
    """Propagate Mean Anomaly from epoch to target JD."""
    P_years = a_AU**1.5
    P_days = P_years * G_DAYS
    n_deg_per_day = 360 / P_days
    d = targetJD - epochJD
    return norm_deg(M0_deg + n_deg_per_day * d)

def heliocentric_from_elements(elem, target_jd):
    """Compute (x, y, z) position (AU) from orbital elements."""
    
    a = elem['a']
    e = elem['e']
    i = deg2rad(elem['i'])
    Omega = deg2rad(elem['Omega'])
    omega = deg2rad(elem['omega'])
    
    Mdeg = propagate_m_deg(elem['M0'], elem['epochJD'], target_jd, a)
    M = deg2rad(Mdeg)
    
    E = solve_kepler(M, e)
    
    nu = 2 * math.atan2(math.sqrt(1+e)*math.sin(E/2), math.sqrt(1-e)*math.cos(E/2))
    r = a * (1 - e * math.cos(E)) # AU

    # Orbital plane coordinates (perihelion-based)
    x_op = r * math.cos(nu)
    y_op = r * math.sin(nu)

    # 1. Rotate by Argument of Perihelion (omega) around z' (orbital Z-axis)
    cosw, sinw = math.cos(omega), math.sin(omega)
    x1 = x_op * cosw - y_op * sinw
    y1 = x_op * sinw + y_op * cosw
    
    # 2. Rotate by Inclination (i) around x' (line of nodes)
    cosi, sini = math.cos(i), math.sin(i)
    x2 = x1
    y2 = y1 * cosi
    z2 = y1 * sini
    
    # 3. Rotate by Longitude of the Ascending Node (Omega) around Z (ecliptic North)
    cosO, sinO = math.cos(Omega), math.sin(Omega)
    x = x2 * cosO - y2 * sinO
    y = x2 * sinO + y2 * cosO
    z = z2
    
    # VPYTHON FIX: Map astronomical (X, Y, Z) to VPython's axes for standard view:
    # Astronomical X (Vernal Equinox) -> VPython X
    # Astronomical Y (90 deg in ecliptic) -> VPython Y
    # Astronomical Z (North Ecliptic Pole) -> VPython Z
    # Since the calculation above yields the correct J2000 (X, Y, Z) coordinates, 
    # we use them directly to form the vector.
    return vector(x, y, z) 

# --- Orbital Elements (J2000-like) ---
ELEMENTS = {
    # Radius (r) is highly exaggerated and based on log(a) for visual scale
    'Mercury': {'a':0.387, 'e':0.205, 'i':7.00, 'Omega':48.33, 'omega':29.12, 'M0':174.79, 'epochJD':J2000_JD, 'color':color.gray(0.7)},
    'Venus':   {'a':0.723, 'e':0.006, 'i':3.39, 'Omega':76.68, 'omega':54.85, 'M0':50.44, 'epochJD':J2000_JD, 'color':color.orange},
    'Earth':   {'a':1.000, 'e':0.016, 'i':0.00, 'Omega':-11.26,'omega':114.20,'M0':357.51, 'epochJD':J2000_JD, 'color':color.green},
    'Mars':    {'a':1.523, 'e':0.093, 'i':1.85, 'Omega':49.57, 'omega':286.46,'M0':19.41, 'epochJD':J2000_JD, 'color':color.red},
    'Jupiter': {'a':5.203, 'e':0.048, 'i':1.30, 'Omega':100.55,'omega':274.19,'M0':19.65, 'epochJD':J2000_JD, 'color':color.yellow},
    'Saturn':  {'a':9.537, 'e':0.054, 'i':2.48, 'Omega':113.71,'omega':338.71,'M0':317.51,'epochJD':J2000_JD, 'color':color.magenta},
    'Uranus':  {'a':19.19, 'e':0.047, 'i':0.76, 'Omega':74.22, 'omega':96.73, 'M0':142.26,'epochJD':J2000_JD, 'color':color.cyan},
    'Neptune': {'a':30.06, 'e':0.008, 'i':1.76, 'Omega':131.72,'omega':273.24,'M0':259.90,'epochJD':J2000_JD, 'color':color.blue},
    'Pluto':   {'a':39.48, 'e':0.248, 'i':17.14, 'Omega':110.29,'omega':113.83,'M0':14.53, 'epochJD':J2000_JD, 'color':color.white}
}

# --- VPython Objects ---

# Sun
sun = sphere(pos=vector(0,0,0), radius=1.5, color=color.yellow, emissive=True)


# Solar System bodies dictionary for VPython objects
planets = {}

# Calculate a visually scaled radius (r_scaled) for the planets
def get_scaled_radius(a_au):
    """Logarithmically scaled radius for visual visibility (0.4 AU to 1.5 AU)."""
    # FIX: Ensure a minimum size for the inner planets
    return 0.4 + math.log10(a_au) * 0.4

# Create the planets and draw their orbits
for name, elem in ELEMENTS.items():
    # Initial position at J2000
    pos_j2000 = heliocentric_from_elements(elem, J2000_JD)
    
    # Create the planet sphere
    planet_radius = get_scaled_radius(elem['a'])
    
    planets[name] = sphere(
        pos=pos_j2000, 
        radius=planet_radius,
        color=elem['color'], 
        make_trail=True,  # Set to True for a live path trail
        trail_type="points",
        interval=50,
        retain=50,
        label=label(text=name, box=False, opacity=0, height=10, xoffset=planet_radius * 1.5)
    )
    # Store the element data in the VPython object
    planets[name].elem = elem
    
    # Create the static orbit (curve)
    planets[name].orbit = curve(color=elem['color'] * 0.4, radius=0.01)
    
# --- Drawing the Orbits (Static Curves) ---

def draw_orbits(elements):
    for name, elem in elements.items():
        orbit_curve = planets[name].orbit
        orbit_curve.clear()
        
        # Sample the orbit using True Anomaly (nu)
        for nu_deg in range(0, 361, 5):
            nu = deg2rad(nu_deg)
            a, e = elem['a'], elem['e']
            
            # Radius in orbital plane
            r = (a * (1 - e*e)) / (1 + e * math.cos(nu))
            
            # Orbital-plane coordinates
            x_op = r * math.cos(nu)
            y_op = r * math.sin(nu)
            
            # Rotation angles
            omega = deg2rad(elem['omega'])
            i = deg2rad(elem['i'])
            Omega = deg2rad(elem['Omega'])

            # Rotation calculations (same as in heliocentric_from_elements)
            x1 = x_op * math.cos(omega) - y_op * math.sin(omega)
            y1 = x_op * math.sin(omega) + y_op * math.cos(omega)
            
            x2 = x1
            y2 = y1 * math.cos(i)
            z2 = y1 * math.sin(i)
            
            x3 = x2 * math.cos(Omega) - y2 * math.sin(Omega)
            y3 = x2 * math.sin(Omega) + y2 * math.cos(Omega)
            z3 = z2
            
            orbit_curve.append(pos=vector(x3, y3, z3))

draw_orbits(ELEMENTS)


# --- Simulation Loop ---

def update_sim(sim_days):
    """Updates the position of all VPython objects based on sim_days."""
    current_jd = J2000_JD + sim_days
    
    for name, planet in planets.items():
        # Get position for current time
        new_pos = heliocentric_from_elements(planet.elem, current_jd)
        
        # Update VPython object position
        planet.pos = new_pos
        planet.label.pos = new_pos
        
    # Update HUD (using scene.caption)
    current_date = datetime.datetime(2000, 1, 1, 12, 0, 0) + datetime.timedelta(days=sim_days)
    
    # Format speed text
    speed_text = f"{time_scale}x FWD" if time_scale > 0 else (f"{-time_scale}x REV" if time_scale < 0 else 'Paused')
    
    scene.caption = f"""
    <h3>Epoch-Aligned Solar System (J2000 Elements)</h3>
    Date: {current_date.strftime('%Y-%m-%d %H:%M:%S')} UTC<br>
    Speed: <b>{speed_text}</b>
    """

def run_sim():
    global sim_days, last_time, time_scale
    
    # Limit frame rate to 30 FPS
    rate(30)
    
    # Calculate time passed since last frame
    now = datetime.datetime.now()
    # FIX: Use .total_seconds() for accurate floating-point duration
    delta_t = (now - last_time).total_seconds() 
    last_time = now
    
    # Advance simulation time
    # delta_t is in seconds. Convert to days by dividing by 86400
    sim_days += (delta_t / 86400.0) * time_scale
    
    # Update all object positions
    update_sim(sim_days)
    
# --- Controls (Widgets) ---

def set_speed(s):
    global time_scale
    time_scale = s.value
    
def go_to_now():
    global sim_days, time_scale
    now = datetime.datetime.now()
    sim_days = days_since_j2000(now)
    time_scale = 1 # Set to real-time after jump
    update_sim(sim_days)

def set_date_button():
    global sim_days, time_scale
    
    # Read text from the input element
    date_str = date_input.text
    try:
        # datetime.strptime can parse the YYYY-MM-DD format
        # We assume UTC noon (12:00:00) for the input day
        target_date = datetime.datetime.strptime(date_str, '%Y-%m-%d').replace(hour=12)
        target_jd = jd_of_date(target_date)
        sim_days = target_jd - J2000_JD
        time_scale = 0 # Pause after jumping to date
    except ValueError:
        scene.caption += "\n<b style='color:red;'>Invalid date format! Use YYYY-MM-DD.</b>"
    
    # Ensure the scene updates immediately after setting the new date
    update_sim(sim_days)


# Create Controls UI
scene.append_to_caption('\n<hr>')
scene.append_to_caption('<h3>Simulation Controls</h3>')

# Speed Slider
scene.append_to_caption('Speed Factor (1x to 500k x): ')
slider(bind=set_speed, min=-500000, max=500000, value=time_scale, step=1000)

# NOW Button
button(text='Set to NOW (1x)', bind=go_to_now)

# Date Input Field and Button
scene.append_to_caption(f"""<br>Go to Date (YYYY-MM-DD): <input type="text" id="dateInput" value="{datetime.datetime.now().strftime('%Y-%m-%d')}">""")
# Get a reference to the HTML element created by the scene.append_to_caption call

button(text='Go to Date (and Pause)', bind=set_date_button)
scene.append_to_caption('\n<hr>')

# --- Kick off the simulation ---

# Set initial time to 'now'
go_to_now()

# Start the VPython animation loop
scene.loop = run_sim

while True: 
    sleep (1)
