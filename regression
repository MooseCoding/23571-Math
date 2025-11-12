import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# === Step 1: Enter your discrete data points ===
# Replace these with your actual data
x = np.array([0.3,0.5,0.5,0.2,0.1,0.65]) # These are your hood positions
y = np.array([950,900,850,850,850,1220]) # These are your flywheel powers
f = np.array([59, 59,49,32,15,107])  # These are your distances

# === Step 2: Choose polynomial degree ===
degree = 2 # Try 2 for quadratic, 3 for cubic, etc.

# === Step 3: Prepare data ===
X = np.column_stack((x, y))

# Create polynomial feature expansion
poly = PolynomialFeatures(degree=degree, include_bias=False)
X_poly = poly.fit_transform(X)

# === Step 4: Fit regression model ===
model = LinearRegression()
model.fit(X_poly, f)

# === Step 5: Display the regression equation ===
terms = poly.get_feature_names_out(['x', 'y'])
coefs = model.coef_
intercept = model.intercept_

print("\nPolynomial Regression Equation:")
equation = f"f(x, y) = {intercept:.4f}"
for term, coef in zip(terms, coefs):
    sign = "+" if coef >= 0 else "-"
    if coef >= 0.00001:
        equation += f" {sign} {abs(coef):.4f}*{term}" 
print(equation)

# === Step 6: Model evaluation ===
f_pred = model.predict(X_poly)
r2 = model.score(X_poly, f)
print(f"R² (Goodness of fit): {r2:.4f}")

# === Step 7: Plot data and regression surface ===
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Scatter plot of data
ax.scatter(x, y, f, color='blue', label='Data Points')

# Create meshgrid for regression surface
x_surf, y_surf = np.meshgrid(np.linspace(min(x), max(x), 40),
                             np.linspace(min(y), max(y), 40))
X_surf = np.column_stack((x_surf.ravel(), y_surf.ravel()))
X_surf_poly = poly.transform(X_surf)
f_surf = model.predict(X_surf_poly).reshape(x_surf.shape)

# Plot the regression surface
ax.plot_surface(x_surf, y_surf, f_surf, alpha=0.5, color='orange')

# Labels
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('f(x, y)')
ax.set_title(f'Polynomial Regression Surface (degree={degree})')
ax.legend()

plt.show()
