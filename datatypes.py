from dataclasses import dataclass
import math
import copy

G = 6.67408e-11  # gravitational constant (you can scale this for visualization)

@dataclass
class OrderedPair:
    """
    A simple 2-D point or vector with named coordinates.

    Used to represent positions, velocities, and accelerations of stars
    in the simulation. Supports access via `.x` and `.y` for readability.
    """
    x: float = 0.0
    y: float = 0.0

    def __iadd__(self, other: "OrderedPair") -> "OrderedPair":
        self.x += other.x
        self.y += other.y
        return self

@dataclass
class Star:
    """
    A celestial body in the simulation.

    Each Star has a position, velocity, acceleration, mass, and radius,
    along with optional RGB color values for visualization.
    The position and motion are expressed as 2D vectors (x, y).
    """
    position: OrderedPair | None = None
    velocity: OrderedPair | None = None
    acceleration: OrderedPair | None = None
    mass: float = 0.0
    radius: float = 0.0
    red: int = 255
    green: int = 255
    blue: int = 255

@dataclass
class Universe:
    """
    A square universe of given width containing a list of stars.

    The universe defines the simulation space. Its width represents the
    side length of a square region with corners at (0, 0) and (width, width).
    """
    width: float = 0.0
    stars: list[Star] = None

    def in_field(self, p: OrderedPair) -> bool:
        """
        Check if a given point is within the bounds of the universe.

        Args:
            p (OrderedPair): The point to check.
        Returns:
            bool: True if the point lies within [0, width] x [0, width],
                  False otherwise.
        """
        return 0 <= p.x <= self.width and 0 <= p.y <= self.width


@dataclass
class Quadrant:
    """
    A square subregion of the universe given by its lower-left corner (x, y) and width.
    """
    x: float = 0.0
    y: float = 0.0
    width: float = 0.0


@dataclass
class Node:
    """
    A quadtree node. Internal nodes store a dummy 'star' for the center of mass and have children.
    Leaf nodes may store a real star (or be empty) and have no children.
    By convention, child quadrants are ordered [NW, NE, SW, SE].
    """
    sector: Quadrant | None = None
    children: list["Node"] | None = None
    star: Star | None = None

    def is_leaf(self) -> bool:
        """
        Check if the node is a leaf (i.e., is child-free).

        Args:
            None
        Returns:
            bool: True if the node has no children, False otherwise.
        """
        return self.children is None or len(self.children) == 0

    def insert(self, s: Star) -> None:
        """
        Insert a star into the quadtree rooted at this node.

        If the node is an empty leaf, the star is stored here directly.
        If the node is an occupied leaf, both the existing and new star are
        pushed down into child nodes. Internal nodes update their center-of-mass
        dummy star and recurse into the appropriate child.

        Args:
            s (Star): The star to insert into the quadtree.
        Returns:
            None
        """
        current_node = self
        # if current_node has no star insert star & exit
        if current_node.star == None:
            current_node.star = s
            return
        
        temp = copy.deepcopy(current_node.star)
        current_node.star.position = center_of_gravity(current_node.star, s)
        current_node.star.mass += s.mass        
        
        # iterate until reach leaf, update nodes along the way
        while not current_node.is_leaf():
            current_node = current_node.find_child(s)
            
            # if star is None then insert star and exit
            if current_node.star == None:
                current_node.star = s
                return
            
            # make copy of star before updating in case star is a leaf and needs to be moved into child
            temp = copy.deepcopy(current_node.star)
            current_node.star.position = center_of_gravity(current_node.star, s)
            current_node.star.mass += s.mass

        # current node has a star, so create children and check if node star and inserted star are  
        #   in different children / quadrants
        current_node.create_children()    
        temp_star_child = current_node.find_child(temp)
        insert_star_child = current_node.find_child(s)
        
        # continue creating children until original star and inserted star are in different quadrants
        while temp_star_child == insert_star_child:
            current_node = temp_star_child
            # set disputed quadrant star to summarize both original star and inserted star
            current_node.star = Star()
            current_node.star.position = center_of_gravity(temp, s)
            current_node.star.mass = temp.mass + s.mass

            # try again to split stars into different quadrants
            current_node.create_children()              
            temp_star_child = current_node.find_child(temp)
            insert_star_child = current_node.find_child(s)
            
        temp_star_child.star = temp 
        insert_star_child.star = s

    def create_children(self) -> None:
        """
        Subdivide this node into four child nodes corresponding to the four
        quadrants of its sector.

        Children are created in the fixed order:
            0: Northwest (NW)
            1: Northeast (NE)
            2: Southwest (SW)
            3: Southeast (SE)
        Each child is an empty leaf node whose sector covers one quarter of
        this node's sector.

        Args:
            None
        Returns:
            None. Modifies self.children in place.
        """
        quad_width = self.sector.width / 2
        NW = Quadrant(self.sector.x, self.sector.y + quad_width, quad_width)
        NE = Quadrant(self.sector.x + quad_width, self.sector.y + quad_width, quad_width)
        SW = Quadrant(self.sector.x, self.sector.y, quad_width)
        SE = Quadrant(self.sector.x + quad_width, self.sector.y, quad_width)
        self.children = [Node(NW), Node(NE), Node(SW), Node(SE)]

    def find_child(self, s: Star) -> 'Node':
        """
        Determine which child node's sector contains the given star.

        Args:
            s (Star): The star whose position determines the target quadrant.
        Returns:
            Node: The child node whose sector would contain star s.
        """
        for child in self.children:
            if child.sector.x <= s.position.x and child.sector.x + child.sector.width > s.position.x and (
                child.sector.y <= s.position.y and child.sector.y + child.sector.width > s.position.y):
                return child
            
        return None

    def calculate_net_force(self, s: Star, theta: float) -> OrderedPair:
        """
        Compute the net gravitational force exerted on a star by all bodies
        in the subtree rooted at this node, using the Barnes-Hut criterion.

        If this node is a sufficiently far-away cluster (sector width / distance
        to center of mass < theta), it is treated as a single body at its
        center of mass. Otherwise, the force contributions of each child are
        summed recursively.

        Args:
            s (Star): The star on which the net force is being computed.
            theta (float): The Barnes-Hut approximation parameter. Larger values
                           trade accuracy for speed; theta=0 gives exact forces.
        Returns:
            OrderedPair: The x- and y- components of the net gravitational force
                         acting on s.
        """
        # return no force if this is a root, it's the same star or stars on on top of each other
        if self.star == None or (self.is_leaf() and self.star is s) or distance(
                self.star.position, s.position) == 0:
            return OrderedPair()
        # if leaf only need to compute with 1 other star
        elif self.is_leaf():
            return compute_force(self.star, s)

        # if distance < theta, treat children as current star
        if distance(self.star.position, s.position) != 0 and self.sector.width / distance(
                self.star.position, s.position) < theta:
            return compute_force(self.star, s)
        # otherwise calculate net force from each child star
        else:
            total_force = OrderedPair()
            for child in self.children:
                total_force += child.calculate_net_force(s, theta)
            return total_force

@dataclass
class QuadTree:
    """
    A wrapper around the root node of a Barnes-Hut quadtree.

    Provides an interface for inserting stars, building the spatial tree,
    and calculating net gravitational forces using hierarchical aggregation.
    """
    root: Node | None = None

    def insert(self, s: Star) -> None:
        """
        Insert a star into the quadtree via the root node.

        Args:
            s (Star): The star to insert.
        Returns:
            None
        """
        self.root.insert(s)

# To prevent circular import issues, we define these functions here.

def center_of_gravity(*stars: Star) -> OrderedPair:
    """
    Compute the center of gravity of an arbitrary number of stars.

    Uses the mass-weighted average of positions:
        x_cg = sum(m_i * x_i) / sum(m_i)
        y_cg = sum(m_i * y_i) / sum(m_i)

    Args:
        *stars (Star): A variadic sequence of Star objects, each with a
                       position (OrderedPair) and a mass (float).
    Returns:
        OrderedPair: The (x, y) coordinates of the collective center of gravity.
    """
    total_px = 0.0
    total_py = 0.0
    count = 0
    
    for star in stars:
        total_px += star.position.x * star.mass
        total_py += star.position.y * star.mass
        count += star.mass
    
    return OrderedPair(round(total_px / count, 3), round(total_py / count, 3))  


def compute_force(s1: Star, s2: Star) -> OrderedPair:
    """
    Compute the gravitational force exerted by s1 on s2.

    Uses Newton's law of universal gravitation: F = G * m1 * m2 / r^2,
    decomposed into x and y components along the displacement vector from
    s2 to s1.

    Args:
        s1 (Star): The star exerting the gravitational force.
        s2 (Star): The star experiencing the gravitational force.
    Returns:
        OrderedPair: The x- and y- components of the force on s2 due to s1.
    """
    d = distance(s1.position, s2.position)
    F = G * s1.mass * s2.mass / (d * d)  
    
    delta = (s1.position.x - s2.position.x, s1.position.y - s2.position.y)
    force = OrderedPair(F * (delta[0] / d), F * (delta[1] / d))
    return force

def distance(p1: OrderedPair, p2: OrderedPair) -> float:
    """
    Compute the Euclidean distance between two points.

    Args:
        p1 (OrderedPair): The first point.
        p2 (OrderedPair): The second point.
    Returns:
        float: The straight-line distance between p1 and p2.
    """
    dx, dy = (p1.x - p2.x, p1.y - p2.y)
    return math.sqrt(dx * dx + dy * dy)
