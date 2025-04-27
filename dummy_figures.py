import math
import csv
import sys
import matplotlib.pyplot as plt
import numpy as np # Import numpy for linspace

def generate_circle_points(radius, num_segments):
    """
    Calculates the coordinates of points evenly distributed on a circle.

    Args:
        radius (float): The radius of the circle.
        num_segments (int): The number of points (segments) to calculate.

    Returns:
        list: A list of [x, y] coordinate pairs.
              Returns an empty list if num_segments is not positive.
    """
    if num_segments <= 0:
        print("Error: Number of segments must be positive.", file=sys.stderr)
        return []
    if radius < 0:
        print("Warning: Radius is negative. Using absolute value.", file=sys.stderr)
        radius = abs(radius)

    points = []
    angle_step = 2 * math.pi / num_segments

    for i in range(num_segments):
        angle = i * angle_step
        # Calculate x and y using trigonometry
        # x = r * cos(theta)
        # y = r * sin(theta)
        x = radius * math.cos(angle)
        y = radius * math.sin(angle)
        points.append([x, y])

    # Optional: Plotting the generated circle points
    # x_coords = [p[0] for p in points]
    # y_coords = [p[1] for p in points]
    # plt.figure()
    # plt.plot(x_coords, y_coords, 'bo-')
    # plt.title(f'Circle R={radius}, N={num_segments}')
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.axis('equal')
    # plt.grid(True)
    # plt.show()

    return points

def generate_parabola_points(num_points):
    """
    Calculates the coordinates of points on the parabola y = x^2
    for x ranging from -1 to 1.

    Args:
        num_points (int): The number of points to calculate.

    Returns:
        list: A list of [x, y] coordinate pairs.
              Returns an empty list if num_points is less than 2.
    """
    if num_points < 2:
        print("Error: Number of points must be at least 2 for a parabola range.", file=sys.stderr)
        return []

    points = []
    # Generate evenly spaced x values from -1 to 1 (inclusive)
    x_values = np.linspace(-1, 1, num_points)

    for x in x_values:
        y = x**2
        points.append([x, y])

    # Optional: Plotting the generated parabola points
    # x_coords = [p[0] for p in points]
    # y_coords = [p[1] for p in points]
    # plt.figure()
    # plt.plot(x_coords, y_coords, 'ro-')
    # plt.title(f'Parabola y=x^2, N={num_points}')
    # plt.xlabel('x')
    # plt.ylabel('y')
    # plt.grid(True)
    # plt.show()

    return points


def write_points_to_csv(points, filename):
    """
    Writes a list of points (x, y coordinates) to a CSV file.

    Args:
        points (list): A list of [x, y] coordinate pairs.
        filename (str): The name of the CSV file to create.

    Returns:
        bool: True if writing was successful, False otherwise.
    """
    if not filename.lower().endswith('.csv'):
        filename += '.csv'
        print(f"Appending .csv extension. Filename is now: {filename}")

    try:
        with open(filename, 'w', newline='') as csvfile:
            writer = csv.writer(csvfile)
            # Write header row
            writer.writerow(['x', 'y'])
            # Write data rows
            writer.writerows(points)
        print(f"Successfully wrote {len(points)} points to '{filename}'")
        return True
    except IOError as e:
        print(f"Error writing to file '{filename}': {e}", file=sys.stderr)
        return False
    except Exception as e:
        print(f"An unexpected error occurred during CSV writing: {e}", file=sys.stderr)
        return False

def main():
    """
    Main function to get user input and orchestrate point generation and CSV writing.
    """
    print("--- Point Generator ---")

    # --- Choose Shape ---
    while True:
        shape_choice = input("Choose shape to generate (c=circle, p=parabola): ").lower()
        if shape_choice in ['c', 'p']:
            break
        else:
            print("Invalid choice. Please enter 'c' or 'p'.")

    generated_points = []
    default_filename = "output.csv" # Default fallback

    # --- Get Parameters and Generate Points ---
    if shape_choice == 'c':
        print("\n--- Circle Parameters ---")
        # Get Radius (R) from user
        while True:
            try:
                r_input = input("Enter the radius (R) of the circle: ")
                radius = float(r_input)
                break # Exit loop if input is a valid float
            except ValueError:
                print("Invalid input. Please enter a numeric value for the radius.")

        # Get Number of Segments (N) from user
        while True:
            try:
                n_input = input("Enter the number of segments/points (N): ")
                num_segments = int(n_input)
                if num_segments <= 0:
                    print("Number of segments must be a positive integer. Please try again.")
                else:
                    break # Exit loop if input is a valid positive integer
            except ValueError:
                print("Invalid input. Please enter an integer value for the number of segments.")

        # Generate points
        generated_points = generate_circle_points(radius, num_segments)
        if generated_points:
            default_filename = f"circle_r{radius}_n{num_segments}.csv"

    elif shape_choice == 'p':
        print("\n--- Parabola Parameters (y=x^2, x from -1 to 1) ---")
        # Get Number of Points (N) from user
        while True:
            try:
                n_input = input("Enter the number of points (N) to generate: ")
                num_points = int(n_input)
                if num_points < 2:
                    print("Number of points must be at least 2. Please try again.")
                else:
                    break # Exit loop if input is a valid integer >= 2
            except ValueError:
                print("Invalid input. Please enter an integer value for the number of points.")

        # Generate points
        generated_points = generate_parabola_points(num_points)
        if generated_points:
            default_filename = f"parabola_n{num_points}.csv"

    # --- Get Filename and Write to CSV ---
    if generated_points: # Proceed only if points were successfully generated
        filename_input = input(f"\nEnter the output CSV filename (press Enter for default: '{default_filename}'): ")
        output_filename = filename_input if filename_input else default_filename
        write_points_to_csv(generated_points, output_filename)
    else:
        print("\nPoint generation failed. No CSV file created.")

if __name__ == "__main__":
    main()
