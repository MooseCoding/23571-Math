import numpy as np
from sklearn.preprocessing import PolynomialFeatures
from sklearn.linear_model import LinearRegression
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

# === Step 1: Enter your discrete data points ===
# Replace these with your actual data
x = np.array([0.3,0.5,0.5,0.2,0.1,0.65]) # These are your hood positions in servo position units 
y = np.array([950,900,850,850,850,1220]) # These are your flywheel powers in the getVelocity() units
f = np.array([59, 59,49,32,15,107])  # These are your distances measured in inches


degree = 2  # Try 2 for quadratic, 3 for cubic, etc.
X = np.column_stack((x, y))
poly = PolynomialFeatures(degree=degree, include_bias=False)
X_poly = poly.fit_transform(X)
model = LinearRegression()
model.fit(X_poly, f)
terms = poly.get_feature_names_out(['x', 'y'])
coefs = model.coef_
intercept = model.intercept_
print("\nPolynomial Regression Equation:")
equation = f"f(x, y) = {intercept:.4f}"
for term, coef in zip(terms, coefs):
    sign = "+" if coef >= 0 else "-"
    equation += f" {sign} {abs(coef):.4f}*{term}"
print(equation)
f_pred = model.predict(X_poly)
r2 = model.score(X_poly, f)
print(f"R² (Goodness of fit): {r2:.4f}")
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(x, y, f, color='blue', label='Data Points')
x_surf, y_surf = np.meshgrid(np.linspace(min(x), max(x), 40),
                             np.linspace(min(y), max(y), 40))
X_surf = np.column_stack((x_surf.ravel(), y_surf.ravel()))
X_surf_poly = poly.transform(X_surf)
f_surf = model.predict(X_surf_poly).reshape(x_surf.shape)
ax.plot_surface(x_surf, y_surf, f_surf, alpha=0.5, color='orange')
ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('f(x, y)')
ax.set_title(f'Polynomial Regression Surface (degree={degree})')
ax.legend()

plt.show()
