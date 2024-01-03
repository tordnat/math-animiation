from manim import *
import numpy as np
import matplotlib.pyplot as plt

def simulate_lorenz_attractor(sigma, rho, beta, time, time_step) -> np.array:
    x, y, z = 1.0, 1.0, 1.0
    timeseries = []

    steps = int(time / time_step)
    for _ in range(steps): # Euler's method
        dx = sigma * (y - x)
        dy = x * (rho - z) - y
        dz = x * y - beta * z

        x = x + dx * time_step
        y = y + dy * time_step
        z = z + dz * time_step

        timeseries.append([x, y, z])

    return np.array(timeseries)

sigma = 10
beta = 8/3
rho = 28
time = 100.0
time_step = 0.001

simulation = simulate_lorenz_attractor(sigma, rho, beta, time, time_step)

class LorenzAttractor(ThreeDScene):
    def construct(self):
        # Setup axes
        axes = ThreeDAxes(
            x_range=[-50, 50, 5],
            y_range=[-50, 50, 5],
            z_range=[-50, 50, 5]
        )

        # Generate points for the Lorenz Attractor
        points = simulate_lorenz_attractor(sigma, rho, beta, time, time_step)
        graph = axes.plot_line_graph(
            x_values=points[:,0], 
            y_values=points[:,1], 
            z_values=points[:,2], 
            add_vertex_dots=False
        )

        # Set the initial camera orientation
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)

        # Add axes to the scene
        self.add(axes)
        self.play(Create(graph), runtime=10)


class Lorenz2Attractor(ThreeDScene):
    def construct(self):
        # Setup axes
        axes = ThreeDAxes(
            x_range=[-50, 50, 5],
            y_range=[-50, 50, 5],
            z_range=[-50, 50, 5]
        )

        # Generate points for the Lorenz Attractor
        points = simulate_lorenz_attractor(sigma, rho, beta, time, time_step)
        graph = axes.plot_line_graph(
            x_values=points[:,0], 
            y_values=points[:,1], 
            z_values=points[:,2], 
            add_vertex_dots=False
        )

        # Set the initial camera orientation
        self.set_camera_orientation(phi=75 * DEGREES, theta= 30*DEGREES)

        # Create a ValueTracker
        tracker = ValueTracker(0)

        # Define an updater function for the graph
        def update_graph(graph):
            # Get the current value of the tracker
            t = tracker.get_value()
            # Determine the portion of the points to be displayed
            portion = points[:int(t * len(points))]
            # Update the graph to display only this portion
            new_graph = axes.plot_line_graph(
                x_values=portion[:,0], 
                y_values=portion[:,1], 
                z_values=portion[:,2], 
                add_vertex_dots=False
            )
            graph.become(new_graph)

        # Apply the updater to the graph
        graph.add_updater(update_graph)

        # Add axes and graph to the scene
        self.add(axes)
        self.add(graph)

        # Animate the ValueTracker over 10 seconds
        self.play(tracker.animate.set_value(1), run_time=10)

        # Remove the updater
        graph.remove_updater(update_graph)

        # [Any other animations or waits]



""" Testing with matplotlib
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Set the same scale for all axes and adjust the viewing angle
ax.set_xlim(-30, 30)
ax.set_ylim(-30, 30)
ax.set_zlim(-30, 30)
ax.plot(simulation[:, 0], simulation[:,1], simulation[:,2])

ax.view_init(elev=30, azim=45)  # Adjust the viewing angle

plt.show()
"""