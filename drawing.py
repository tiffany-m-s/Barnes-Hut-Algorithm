import pygame
from datatypes import Universe
from engine import barnes_hut


def animate_system(
    time_points: list[Universe],
    canvas_width: int,
    frequency: int,
    scaling_factor: float
) -> list[pygame.Surface]:
    """
    Convert a sequence of Universe snapshots into a list of Pygame surfaces.

    Iterates over time_points and renders every frequency-th universe to a
    canvas, collecting the resulting surfaces for downstream use (e.g., saving
    as a video).

    Args:
        time_points (list[Universe]): Ordered list of Universe states produced
                                      by the Barnes-Hut simulation.
        canvas_width (int): The side length in pixels of the square canvas.
        frequency (int): Render one frame for every this many time steps.
                         For example, frequency=100 renders frames at indices
                         0, 100, 200, ...
        scaling_factor (float): Multiplier applied to star radii when drawing,
                                so stars remain visible relative to the canvas.
    Returns:
        list[pygame.Surface]: A list of rendered Pygame surfaces, one per
                              sampled time step, in chronological order.
    """
    images = []

    for i, universe in enumerate(time_points):
        if i % frequency == 0:
            img = draw_to_canvas(universe, canvas_width, scaling_factor)
            images.append(img)

    return images


def draw_to_canvas(
    universe: Universe,
    canvas_width: int,
    scaling_factor: float
) -> pygame.Surface:
    """
    Render a single Universe snapshot onto a Pygame surface and return it.

    Each star is drawn as a filled circle whose screen position is scaled
    proportionally to the universe width, and whose radius is scaled by
    scaling_factor. Stars with no position data are skipped.

    Args:
        universe (Universe): The universe snapshot to render, containing a
                             width and a list of Star objects.
        canvas_width (int): The side length in pixels of the square output surface.
        scaling_factor (float): Multiplier applied to each star's radius before
                                mapping it to pixel space, so stars remain
                                visible at simulation scales.
    Returns:
        pygame.Surface: A canvas_width x canvas_width surface with a black
                        background and all stars drawn in their respective colors.
    """
    surface = pygame.Surface((canvas_width, canvas_width))
    surface.fill((0, 0, 0))

    if universe.stars is None:
        return surface

    for star in universe.stars:
        color = (star.red, star.green, star.blue)
        cx = int((star.position.x / universe.width) * canvas_width)
        cy = int((star.position.y / universe.width) * canvas_width)
        r = max(1, int(scaling_factor * (star.radius / universe.width) * canvas_width))
        pygame.draw.circle(surface, color, (cx, cy), r)

    return surface
