import numpy as np
import math
import os

# ---------- Configuration ----------
class_name = "Aimbot"
output_path = f"{class_name}.java"   # writes to current directory

# numeric scan parameters
k_start = 6.0
k_end = 204.0
k_step = 0.1

y_min = 750.0
y_max = 2000.0

x_min = 0.0
x_max = 0.7 

# ---------- Solver ----------
def solve_for_x(k, y):
    """
    Solve the quadratic for x given k and y.
    Returns the in-range root closest to k, or None if no real root exists.
    """
    a = -283.3333
    b = 98.7865 + 0.1838 * y
    c = -200.1912 + 0.3383 * y - 0.0001 * y**2 - k

    disc = b*b - 4*a*c
    if disc < 0:
        return None

    sqrt_disc = math.sqrt(disc)
    roots = [(-b + sqrt_disc)/(2*a), (-b - sqrt_disc)/(2*a)]

    # Filter roots in range
    in_range = [x for x in roots if x_min <= x <= x_max and math.isfinite(x)]
    if not in_range:
        return None

    # Pick the root closest to desired k
    best_x = min(in_range, key=lambda x: abs((-283.3333*x**2 + (98.7865+0.1838*y)*x + 0.3383*y - 200.1912 - 0.0001*y**2) - k))
    return best_x

# ---------- Formatter ----------
def fmt_val(v):
    return f"{v:.3f}" if np.isfinite(v) else "Double.NaN"

# ---------- Generate lookup points ----------
k_values = np.arange(k_start, k_end + 1e-9, k_step)
y_scan = np.linspace(y_min, y_max, 2000)  # finer scan helps avoid NaNs

points = []
for k in k_values:
    best_pair = (float('nan'), float('nan'))
    best_res = float('inf')
    for y in y_scan:
        x = solve_for_x(k, y)
        if x is not None:
            # Residual for best selection
            res = abs((-283.3333*x**2 + (98.7865+0.1838*y)*x + 0.3383*y - 200.1912 - 0.0001*y**2) - k)
            if res < best_res:
                best_res = res
                best_pair = (x, y)
    points.append((round(best_pair[0],3), round(best_pair[1],3)))
    print(f"Hood: {best_pair[0]:.3f}, Flywheel: {best_pair[1]:.3f}, Distance: {k:.3f}")

# ---------- Build Java Source ----------
java_lines = [
    "/**",
    " * Auto-generated Aimbot lookup table.",
    " * Each entry: { hood, flywheel } rounded to 3 decimals.",
    " */",
    f"public class {class_name} " + "{",
    "    public static final double[][] points = {"
]

for x, y in points:
    java_lines.append(f"        {{ {fmt_val(x)}, {fmt_val(y)} }},")
java_lines.append("    };")
java_lines.extend([
    "",
    "    public static double[] lerp(double t, double[] low, double[] high) {",
    "        double hood = low[0] + t * (high[0] - low[0]);",
    "        double flywheel = low[1] + t * (high[1] - low[1]);",
    "        return new double[] { hood, flywheel };",
    "    }",
    "",
    f"    public static double[] getValues(double distance) {{",
    f"        double k_step = {k_step};",
    "        if (Double.isNaN(distance) || k_step <= 0.0) return new double[] { Double.NaN, Double.NaN };",
    "        int n = points.length;",
    "        if (n == 0) return new double[] { Double.NaN, Double.NaN };",
    "        int low = (int) Math.floor(distance / k_step);",
    "        int high = (int) Math.ceil(distance / k_step);",
    "        low = Math.max(0, Math.min(low, n - 1));",
    "        high = Math.max(0, Math.min(high, n - 1));",
    "        if (high == low) return points[low];",
    "        double t = (distance - low * k_step) / k_step;",
    "        return lerp(t, points[low], points[high]);",
    "    }",
    "",
    f"    private {class_name}() {{}}",
    "}"
])

java_source = "\n".join(java_lines)

# ---------- Write Java file ----------
with open(output_path, "w", newline="\n") as f:
    f.write(java_source)

print(f"Wrote Java class to: {os.path.abspath(output_path)}")
print(f"Rows written: {len(points)}")
