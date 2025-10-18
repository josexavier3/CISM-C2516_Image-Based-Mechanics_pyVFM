import numpy as np

def compute_Q_matrix(E1, E2, nu12, G12):
    """
    Compute the reduced stiffness matrix Q for plane stress conditions
    
    Parameters:
    E1, E2: Young's moduli in principal directions (GPa)
    nu12: Poisson's ratio 
    G12: Shear modulus (GPa)
    
    Returns:
    Q: 3x3 reduced stiffness matrix (GPa)
    """
    
    # Calculate nu21 using symmetry condition
    nu21 = nu12 * E2 / E1
    
    # Calculate the determinant term
    delta = 1 - nu12 * nu21
    
    # Compute Q matrix components
    Q11 = E1 / delta
    Q22 = E2 / delta  
    Q12 = nu12 * E2 / delta  # or equivalently: nu21 * E1 / delta
    Q66 = G12
    
    # Assemble the Q matrix
    Q = np.array([
        [Q11, Q12, 0],
        [Q12, Q22, 0],
        [0,   0,   Q66]
    ])
    
    return Q, nu21, delta

# Material properties from your data
print("========================")
E1 = 15.1  # GPa
E2 = 1.91  # GPa  
nu12 = 0.471
G12 = 1.109  # GPa

print(f"Input parameters:")
print(f"E1 = {E1} GPa")
print(f"E2 = {E2} GPa") 
print(f"ν12 = {nu12}")
print(f"G12 = {G12} GPa")

Q_actual, nu21_actual, delta_actual = compute_Q_matrix(E1, E2, nu12, G12)

print(f"\nCalculated values:")
print(f"ν21 = {nu21_actual:.4f}")
print(f"δ = 1 - ν12×ν21 = {delta_actual:.4f}")

print(f"\nComputed Q matrix (GPa):")
print(f"Q = [{Q_actual[0,0]:.3f}  {Q_actual[0,1]:.3f}   {Q_actual[0,2]:.3f}]")
print(f"    [{Q_actual[1,0]:.3f}  {Q_actual[1,1]:.3f}   {Q_actual[1,2]:.3f}]")
print(f"    [{Q_actual[2,0]:.3f}   {Q_actual[2,1]:.3f}   {Q_actual[2,2]:.3f}]")
