import sys
import argparse
import numpy as np
import os
import json

def main():
    parser = argparse.ArgumentParser(description="Predict DBBA monomers and J (extrapolated MFH) for a 7-AGNR.")
    parser.add_argument("length_nm", type=float, help="Experimental length of the nanoribbon (in nm)")
    args = parser.parse_args()

    input_length = args.length_nm
    print(f"\n--- Analyzing experimental length: {input_length:.3f} nm ---")

    try:
        import sisl
        from sisl_hubbard_7agnr import build_7agnr_geometry, calculate_J_mfh
    except ImportError:
        print("Error: Could not import 'sisl' or 'sisl_hubbard_7agnr'. Run this script in your correct conda/venv environment.")
        sys.exit(1)

    # 1. Fast determination of closest L
    print("Determining geometry lengths...")
    L_to_length = {}
    for L in range(2, 12):
        geom = build_7agnr_geometry(L)
        x_coords = geom.xyz[:, 0]
        L_to_length[L] = (x_coords.max() - x_coords.min()) / 10.0  # Convert to nm
        
    closest_L = min(L_to_length.keys(), key=lambda L: abs(L_to_length[L] - input_length))
    
    # We will compute MFH only up to the closest L (with a minimum of L=3 to allow exponential fitting)
    max_L_to_compute = max(closest_L, 3)
    
    print(f"Target length matches closely with L={closest_L}.")
    print(f"Calculating MFH data up to L={max_L_to_compute} to perform the extrapolation...")

    lengths_nm = []
    J_mfh_meV = []
    L_vals_computed = []

    for L in range(2, max_L_to_compute + 1):
        geom = build_7agnr_geometry(L)
        length_nm = L_to_length[L]
        
        try:
            rm = calculate_J_mfh(geom, U_val=3.0)
            J_val = rm['J'] * 1000  # Convert to meV
            
            lengths_nm.append(length_nm)
            J_mfh_meV.append(J_val)
            L_vals_computed.append(L)
            
            print(f"  Computed L={L}: Length = {length_nm:.3f} nm, J_MFH = {J_val:.4e} meV")
        except Exception as e:
            print(f"  Error computing L={L}: {e}")

    lengths_nm = np.array(lengths_nm)
    J_mfh_meV = np.array(J_mfh_meV)
    L_vals_computed = np.array(L_vals_computed)

    # Filter out invalid values
    valid_idx = np.isfinite(J_mfh_meV) & (J_mfh_meV > 0)
    
    if np.sum(valid_idx) < 2:
        print("Error: Not enough valid MFH data points to extrapolate.")
        sys.exit(1)

    # 2. Extrapolate equivalent monomers (continuous L)
    # L scales linearly with length
    coeffs_L = np.polyfit(lengths_nm, L_vals_computed, 1)
    L_extrapolated = np.polyval(coeffs_L, input_length)

    # 3. Extrapolate J_MFH
    # Exponential fit: J = A * exp(-length / xi) => ln(J) = ln(A) - length / xi
    log_J = np.log(J_mfh_meV[valid_idx])
    coeffs_J = np.polyfit(lengths_nm[valid_idx], log_J, 1)
    
    # Calculate expected J for the exact experimental length
    log_J_expected = np.polyval(coeffs_J, input_length)
    J_expected = np.exp(log_J_expected)

    # Calculate J for exactly the closest integer L
    log_J_closest_L = np.polyval(coeffs_J, L_to_length[closest_L])
    J_closest_L = np.exp(log_J_closest_L)

    print("\n================ RESULTS ================")
    print(f"Equivalent monomers (continuous) : {L_extrapolated:.2f} DBBA")
    print(f"Closest monomers (integer L)     : {closest_L} DBBA")
    print("-" * 43)
    print(f"Extrapolated J for exact length of {input_length:.3f} nm:")
    print(f"  -> {J_expected:.4e} meV")
    print("-" * 43)
    print(f"Extrapolated J for exactly L = {closest_L}:")
    print(f"  -> {J_closest_L:.4e} meV")
    print("=========================================")

if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    main()
