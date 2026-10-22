import math
from datatypes import OrderedPair, Universe, QuadTree, Node, Quadrant, Star, distance, compute_force, center_of_gravity
from copy import deepcopy


def barnes_hut(initial_universe: Universe, num_gens: int, time: float, theta: float) -> list[Universe]:
    """
    Run a Barnes-Hut N-body gravitational simulation for a given number of steps.

    Starting from initial_universe, repeatedly applies update_universe() to
    advance the simulation by one time step per generation.

    Args:
        initial_universe (Universe): The starting state of the simulation,
                                     containing all stars and the universe width.
        num_gens (int): The number of simulation steps to run.
        time (float): The duration of each time step in seconds.
        theta (float): The Barnes-Hut approximation parameter passed to
                       calculate_net_force(). Larger values are faster but
                       less accurate; theta=0 gives exact pairwise forces.
    Returns:
        list[Universe]: A list of num_gens+1 Universe objects, where index 0
                        is initial_universe and index i is the state after
                        the i-th simulation step.
    """
    gens = [initial_universe]
    for i in range(num_gens):
        gens.append(update_universe(gens[-1], time, theta))
        
    return gens


def update_universe(current_universe: Universe, time: float, theta: float) -> Universe:
    """
    Advance the universe by one time step using the Barnes-Hut algorithm.

    Builds a quadtree from the current universe, then for each star computes
    its updated acceleration, velocity, and position using that tree and the
    trapezoidal integration helpers.

    Args:
        current_universe (Universe): The current state of the simulation.
        time (float): The duration of the time step in seconds.
        theta (float): The Barnes-Hut approximation parameter. Controls the
                       accuracy vs. speed tradeoff in force calculations.
    Returns:
        Universe: A new Universe object representing the state of the system
                  after one time step.
    """
    new_universe = copy_universe(current_universe)
    tree = generate_quadtree(new_universe)

    # update each star in new_universe using current_universe to prevent affecting fetch
    for i in range(len(current_universe.stars)):
        new_universe.stars[i].acceleration = update_acceleration(new_universe.stars[i], tree, theta)
        new_universe.stars[i].velocity = update_velocity(new_universe.stars[i], time, 
                                                         current_universe.stars[i].acceleration)
        new_universe.stars[i].position = update_position(new_universe.stars[i], time, 
                current_universe.stars[i].acceleration, current_universe.stars[i].velocity)
        
    return new_universe

def generate_quadtree(universe: Universe) -> QuadTree:
    """
    Build a Barnes-Hut quadtree from all stars in the given universe.

    Creates a root node covering the full universe sector, then inserts
    every star into the tree. The resulting tree can be used to efficiently
    approximate gravitational forces via calculate_net_force().

    Args:
        universe (Universe): The current universe, containing a list of stars
                             and a width defining the simulation boundaries.
    Returns:
        QuadTree: A QuadTree whose root covers the full universe and contains
                  all stars from universe.stars.
    """
    tree = QuadTree()
    temp = Quadrant()
    temp.width = universe.width
    tree.root = Node(temp)
    # check if each star is in_field before inserting a copy (to prevent changing variables from
    #   affecting fetch)
    for star in universe.stars:
        if universe.in_field(star.position):
            tree.insert(deepcopy(star))
            
    return tree

G = 6.67408e-11  # gravitational constant (you can scale this for visualization)

# Here are some helper functions for updating star properties.
# They use the quadtree to compute gravitational forces efficiently.
# They are mostly taken from our previous n-body simulation code
# (note that update_acceleration() requires the quadtree).

def update_acceleration(star: Star, q: QuadTree, theta: float) -> OrderedPair:
    """
    Compute the gravitational acceleration acting on a star using the quadtree.

    Queries the quadtree for the net force on the star, then divides by the
    star's mass (F = ma => a = F/m).

    Args:
        star (Star): The star whose acceleration is being updated.
        q (QuadTree): The Barnes-Hut quadtree built from the current universe.
        theta (float): The Barnes-Hut approximation parameter passed through
                       to calculate_net_force().
    Returns:
        OrderedPair: The x- and y- components of the star's new acceleration.
    """
    if q.root is None or star.mass == 0.0:
        return OrderedPair(0.0, 0.0)

    force = (q.root).calculate_net_force(star, theta)

    if star.mass == 0.0:
        return OrderedPair(0.0, 0.0)
    else:
        ax = force.x / star.mass
        ay = force.y / star.mass
        return OrderedPair(ax, ay)

def update_velocity(star: Star, time: float, old_accel: OrderedPair) -> OrderedPair:
    """
    Compute a star's new velocity using the trapezoidal integration rule.

    Averages the star's previous acceleration (old_accel) and its current
    acceleration (star.acceleration) over the time step.

    Args:
        star (Star): The star whose velocity is being updated. Must have a
                     valid velocity and acceleration set.
        time (float): The duration of the time step in seconds.
        old_accel (OrderedPair): The star's acceleration at the start of the
                                 time step (before the update).
    Returns:
        OrderedPair: The x- and y- components of the star's updated velocity.
    """
    star.acceleration = star.acceleration or OrderedPair(0.0, 0.0)
    new_velocity = OrderedPair(0.0, 0.0)

    dvx = 0.5 * (old_accel.x + star.acceleration.x) * time
    dvy = 0.5 * (old_accel.y + star.acceleration.y) * time
    new_velocity = OrderedPair(star.velocity.x + dvx, star.velocity.y + dvy)

    return new_velocity

def update_position(star: Star, time: float, old_accel: OrderedPair, old_vel: OrderedPair) -> OrderedPair:
    """
    Compute a star's new position using constant-acceleration kinematics.

    Applies the equation: delta_x = v0*t + 0.5*a*t^2, using the star's
    velocity and acceleration at the start of the time step.

    Args:
        star (Star): The star whose position is being updated.
        time (float): The duration of the time step in seconds.
        old_accel (OrderedPair): The star's acceleration at the start of the
                                 time step.
        old_vel (OrderedPair): The star's velocity at the start of the time step.
    Returns:
        OrderedPair: The x- and y- coordinates of the star's updated position.
    """
    new_position = OrderedPair(0.0, 0.0)

    dx = 0.5 * old_accel.x * time * time + old_vel.x * time
    dy = 0.5 * old_accel.y * time * time + old_vel.y * time
    new_position = OrderedPair(star.position.x + dx, star.position.y + dy)

    return new_position

def copy_universe(current_universe: Universe) -> Universe:
    """
    Create a shallow-independent copy of a Universe with duplicated Star objects.

    Each Star in the new Universe gets freshly allocated position, velocity,
    and acceleration OrderedPairs, so mutations to the copy will not affect
    the original. Faster than deepcopy for this known structure.

    Args:
        current_universe (Universe): The universe to copy.
    Returns:
        Universe: A new Universe with the same width and independent Star copies.
    """
    new_stars: list[Star] = []

    for s in current_universe.stars:
        new_star = Star(
            position=OrderedPair(s.position.x, s.position.y),
            velocity=OrderedPair(s.velocity.x, s.velocity.y),
            acceleration=OrderedPair(s.acceleration.x, s.acceleration.y),
            mass=s.mass,
            radius=s.radius,
            red=s.red,
            green=s.green,
            blue=s.blue,
        )
        new_stars.append(new_star)

    return Universe(width=current_universe.width, stars=new_stars)
