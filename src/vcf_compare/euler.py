import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.axes import Axes

# Regions
# 1  old_all only          (not new_all, not passing)
# 2  new_all only          (not old_all, not passing)
# 3  old_all & new_all     (shared, failing in both)
# 4  old_pass only         (in old_all but not new_all)
# 5  new_pass only         (in new_all but not old_all)
# 6  old_pass & new_all    (regression — was passing, now failing)
# 7  new_pass & old_all    (improvement — was failing, now passing)
# 8  old_pass & new_pass   (passing in both runs)

# TODO - simple-venn parity
# - set_labels
# - set_colours
# - set font sizes


def plot_pass_fail_euler_diagram(sizes: list[int], ax: Axes | None = None) -> Axes:
    """Plot nested venn diagram / euler diagram showing the intersections of 2 sets and their pass subsets"""
    if len(sizes) != 8:
        raise ValueError(f"Invalid number of sizes {len(sizes)} for pass/fail plot, must be 8!")

    # Geometry:
    # Two large circles (all), two smaller circles (pass) contained within them.
    # pass circles are offset downward so they overlap each other inside the
    # shared region of the two all-circles.

    all_radius = 2.0
    pass_radius = 1.1

    old_all_pos = (-1.2, 0.3)
    new_all_pos = (1.2, 0.3)
    old_pass_pos = (-0.7, -0.4)
    new_pass_pos = (0.7, -0.4)

    label_pos = [
        (-2.2, 1.10),  # inside old_all, outside new_all and both pass circles
        (2.2, 1.10),  # symmetric
        (0.0, 1.10),  # inside both all-circles, above both pass circles
        (-1.2, -0.50),  # inside old_pass, outside new_all
        (1.2, -0.50),  # symmetric
        (-0.5, 0.30),  # inside old_pass & new_all, outside new_pass
        (0.5, 0.30),  # symmetric
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
        (old_all_pos, all_radius, "#7B6FCF"),  # old_all — purple
        (new_all_pos, all_radius, "#2A9D7A"),  # new_all — teal
        (old_pass_pos, pass_radius, "#534AB7"),  # old_pass — dark purple
        (new_pass_pos, pass_radius, "#0F6E56"),  # new_pass — dark teal
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
        (-2.4, old_all_pos[1] + all_radius + 0.08, "old_all", 13, "bold", "normal", "#3C3489"),
        (2.4, new_all_pos[1] + all_radius + 0.08, "new_all", 13, "bold", "normal", "#085041"),
        (
            old_pass_pos[0] - 0.5,
            old_pass_pos[1] - pass_radius - 0.15,
            "old_pass",
            11,
            "normal",
            "italic",
            "#534AB7",
        ),
        (
            new_pass_pos[0] + 0.5,
            new_pass_pos[1] - pass_radius - 0.15,
            "new_pass",
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
