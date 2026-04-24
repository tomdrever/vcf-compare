import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.axes import Axes

# Regions
# 1  a_all only          (not b_all, not passing)
# 2  b_all only          (not a_all, not passing)
# 3  a_all & b_all     (shared, failing in both)
# 4  a_pass only         (in a_all but not b_all)
# 5  b_pass only         (in b_all but not a_all)
# 6  a_pass & b_all    (regression — was passing, now failing)
# 7  b_pass & a_all    (improvement — was failing, now passing)
# 8  a_pass & b_pass   (passing in both runs)

# TODO - simple-venn parity
# - set_labels
# - set_colours
# - set font sizes


def plot_pass_fail_euler_diagram(sizes: list[int], a_name: str, b_name: str, ax: Axes | None = None) -> Axes:
    """Plot nested venn diagram / euler diagram showing the intersections of 2 sets and their pass subsets"""
    if len(sizes) != 8:
        raise ValueError(f"Invalid number of sizes {len(sizes)} for pass/fail plot, must be 8!")

    # Geometry:
    # Two large circles (all), two smaller circles (pass) contained within them.
    # pass circles are offset downward so they overlap each other inside the
    # shared region of the two all-circles.

    all_radius = 2.0
    pass_radius = 1.1

    a_all_pos = (-1.2, 0.3)
    b_all_pos = (1.2, 0.3)
    a_pass_pos = (-0.7, -0.4)
    b_pass_pos = (0.7, -0.4)

    a_all_str = f"{a_name}_all"
    b_all_str = f"{b_name}_all"
    a_pass_str = f"{a_name}_pass"
    b_pass_str = f"{b_name}_pass"

    label_pos = [
        (-2.2, 1.10),  # inside a_all, outside b_all and both pass circles
        (2.2, 1.10),   # symmetric
        (0.0, 1.10),   # inside both all-circles, above both pass circles
        (-1.2, -0.50), # inside a_pass, outside b_all
        (1.2, -0.50),  # symmetric
        (-0.5, 0.30),  # inside a_pass & b_all, outside b_pass
        (0.5, 0.30),   # symmetric
        (0.0, -0.50),  # inside both pass circles
    ]

    # Plot
    if ax is None:
        _, ax = plt.subplots(figsize=(12, 8))

    ax.set_aspect("equal")
    ax.set_xlim(-3.6, 3.6)
    ax.set_ylim(-2.0, 2.4)
    ax.axis("off")

    # Draw circles in order - all-circles first, then pass-circles on top
    circle_specs = [
        (a_all_pos, all_radius, "#7B6FCF"),   # a_all — purple
        (b_all_pos, all_radius, "#2A9D7A"),   # b_all — teal
        (a_pass_pos, pass_radius, "#534AB7"), # a_pass — dark purple
        (b_pass_pos, pass_radius, "#0F6E56"), # b_pass — dark teal
    ]

    for (cx, cy), radius, colour in circle_specs:
        ax.add_patch(
            mpatches.Circle(
                (cx, cy),
                radius,
                facecolor=colour,
                alpha=0.18,
                edgecolor=colour,
                linewidth=2.0,
                zorder=1,
            )
        )

    # Circle name labels
    label_specs = [
        (-2.4, a_all_pos[1] + all_radius + 0.08, a_all_str, 13, "bold", "normal", "#3C3489"),
        (2.4, b_all_pos[1] + all_radius + 0.08, b_all_str, 13, "bold", "normal", "#085041"),
        (
            a_pass_pos[0] - 0.5,
            a_pass_pos[1] - pass_radius - 0.15,
            a_pass_str,
            11,
            "normal",
            "italic",
            "#534AB7",
        ),
        (
            b_pass_pos[0] + 0.5,
            b_pass_pos[1] - pass_radius - 0.15,
            b_pass_str,
            11,
            "normal",
            "italic",
            "#0F6E56",
        ),
    ]

    for x, y, text, size, weight, style, colour in label_specs:
        ax.text(
            x,
            y,
            text,
            ha="center",
            va="bottom",
            fontsize=size,
            fontweight=weight,
            fontstyle=style,
            color=colour,
        )

    # Region count labels
    for region, (px, py) in enumerate(label_pos):
        text = f"{sizes[region]}"
        ax.text(
            px,
            py,
            text,
            ha="center",
            va="center",
            fontsize=12,
            color="#111111",
            bbox=dict(
                boxstyle="round,pad=0.35",
                facecolor="white",
                alpha=0.80,
                edgecolor="none",
            ),
            zorder=2,
        )

    ax.set_title("Old vs new — all / pass breakdown", fontsize=14, pad=30)
    return ax
