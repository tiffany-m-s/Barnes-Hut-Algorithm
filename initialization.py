import math
import random
from datatypes import OrderedPair, Star, Universe  # Adjust this path if needed

G = 6.67408e-11  # gravitational constant
M = 8e36         # mass of black hole

# A Galaxy is just a list of Star objects
Galaxy = list[Star]


def initialize_universe(galaxies: list[Galaxy], w: float) -> Universe:
    """
    Construct an initial Universe from a list of galaxies.

    Flattens all stars from every galaxy into a single list and wraps them
    in a Universe object with the given width.

    Args:
        galaxies (list[Galaxy]): A list of galaxies, where each Galaxy is
                                 itself a list of Star objects.
        w (float): The side length of the square universe in meters.
    Returns:
        Universe: A Universe of width w containing every star from all
                  provided galaxies.
    """
    stars = []
    for galaxy in galaxies:
        for star in galaxy:
            stars.append(star)
    return Universe(width=w, stars=stars)


def initialize_galaxy(num_of_stars: int, r: float, x: float, y: float) -> Galaxy:
    """
    Generate a galaxy of stars orbiting a central black hole.

    Each star is placed at a random distance in [r/2, r] from the center
    (x, y) at a random angle, and given a tangential orbital velocity
    approximated by circular-orbit mechanics. A stationary black hole of
    mass M is appended at the galaxy center.

    Star properties (per star):
        - mass:   1.989e30 kg  (one solar mass)
        - radius: 696,340,000 m (one solar radius)
        - color:  white (255, 255, 255)

    Black hole properties:
        - mass:   M = 8e36 kg
        - radius: 6,963,400,000 m (10x solar radius)
        - color:  blue (0, 0, 255)

    Args:
        num_of_stars (int): The number of orbiting stars to generate,
                            not counting the central black hole.
        r (float): The maximum orbital radius from the galaxy center in meters.
                   Stars are placed at distances uniformly sampled from
                   [r/2, r].
        x (float): The x-coordinate of the galaxy center in meters.
        y (float): The y-coordinate of the galaxy center in meters.
    Returns:
        Galaxy: A list of num_of_stars Star objects plus one black hole Star,
                for a total of num_of_stars + 1 entries.
    """
    galaxy = []

    for _ in range(num_of_stars):
        # First choose distance to center of galaxy and multiply by factor of r
        dist = (random.random() + 1.0) / 2.0 * r

        # Next choose the angle in radians to represent the rotation
        angle = random.random() * 2 * math.pi

        # convert polar coordinates to Cartesian
        pos_x = x + dist * math.cos(angle)
        pos_y = y + dist * math.sin(angle)
        position = OrderedPair(pos_x, pos_y)

        # orbital velocity approximation
        speed = math.sqrt(G * M / dist)
        vel_x = speed * math.cos(angle + math.pi / 2)
        vel_y = speed * math.sin(angle + math.pi / 2)
        velocity = OrderedPair(vel_x, vel_y)

        # create star
        star = Star(
            position=position,
            velocity=velocity,
            acceleration=OrderedPair(0, 0),
            mass=1.989e30,       # set the mass = mass of sun by default
            radius=696340000,    # set the radius equal to radius of sun in m
            red=255,
            green=255,
            blue=255
        )

        galaxy.append(star)

    # Add black hole to the center of the galaxy
    blackhole = Star(
        position=OrderedPair(x, y),
        velocity=OrderedPair(0, 0),
        acceleration=OrderedPair(0, 0),
        mass=M,
        radius=6963400000,
        red=0,
        green=0,
        blue=255
    )

    galaxy.append(blackhole)
    return galaxy


def push(galaxy: Galaxy, v: OrderedPair) -> None:
    """
    Apply a uniform velocity offset to every star in a galaxy.

    Adds the vector v to each star's current velocity, effectively giving
    the entire galaxy a bulk translational motion. Useful for setting up
    colliding galaxies with non-zero initial velocities.

    Args:
        galaxy (Galaxy): The list of Star objects to update.
        v (OrderedPair): The velocity vector (in m/s) to add to each star.
    Returns:
        None. Modifies each star's velocity in place.
    """
    for star in galaxy:
        star.velocity = OrderedPair(star.velocity.x + v.x, star.velocity.y + v.y)