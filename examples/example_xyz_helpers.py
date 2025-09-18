#!/usr/bin/env python3
"""Example usage of XYZ helper methods for third-party library integration."""

import asyncio
import os
from atct import get_species, search_species

async def main():
    """Demonstrate XYZ helper functionality with real API data."""
    print("=== ATcT XYZ Helper Examples ===\n")
    
    # Get species with XYZ data from the API
    print("1. Fetching species with XYZ data from API...")
    try:
        # Try to get a species that likely has XYZ data
        species = await get_species("67-56-1*0", expand_xyz=True)  # Methanol (gas)
        print(f"   Species: {species.name} ({species.formula})")
        print(f"   ATcT ID: {species.atct_id}")
        print(f"   SMILES: {species.smiles}")
    except Exception as e:
        print(f"   Error fetching species: {e}")
        print("   Trying alternative species...")
        # Try a different species if the first one fails
        try:
            species = await get_species("74-82-8*0", expand_xyz=True)  # Methane (gas)
            print(f"   Species: {species.name} ({species.formula})")
            print(f"   ATcT ID: {species.atct_id}")
            print(f"   SMILES: {species.smiles}")
        except Exception as e2:
            print(f"   Error fetching alternative species: {e2}")
            return
    print()
    
    # Check if XYZ data is available
    print("2. Checking XYZ data availability...")
    xyz_data = species.get_xyz_data()
    if xyz_data:
        print(f"   ✓ XYZ data available: {len(xyz_data)} atoms")
        print(f"   Atoms: {xyz_data.atoms}")
        print(f"   Comment: {xyz_data.comment}")
    else:
        print("   ✗ No XYZ data available for this species")
        print("   This might be because:")
        print("   - The species doesn't have XYZ coordinates in the database")
        print("   - The XYZ data is stored as a filename (not expanded)")
        print("   - There was an error parsing the XYZ data")
        return
    print()
    
    # Show XYZ data
    print("3. XYZ coordinate data:")
    xyz_string = species.get_xyz_string()
    print("   " + "\n   ".join(xyz_string.split('\n')))
    print()
    
    # Demonstrate third-party library integration
    print("4. Third-party library integration:")
    
    # Check which libraries are available
    ase_available = False
    rdkit_available = False
    pymatgen_available = False
    
    try:
        import ase
        ase_available = True
        print("   ✓ ASE is available")
    except ImportError:
        print("   ✗ ASE not available (install with: pip install ase)")
    
    try:
        import rdkit
        rdkit_available = True
        print("   ✓ RDKit is available")
    except ImportError:
        print("   ✗ RDKit not available (install with: pip install rdkit)")
    
    try:
        import pymatgen
        pymatgen_available = True
        print("   ✓ pymatgen is available")
    except ImportError:
        print("   ✗ pymatgen not available (install with: pip install pymatgen)")
    print()
    
    # ASE integration
    print("5. ASE (Atomic Simulation Environment) integration:")
    if ase_available:
        try:
            ase_atoms = species.to_ase_atoms()
            print(f"   ✓ Created ASE Atoms object")
            print(f"   Chemical formula: {ase_atoms.get_chemical_formula()}")
            print(f"   Number of atoms: {len(ase_atoms)}")
            print(f"   Positions shape: {ase_atoms.get_positions().shape}")
            print(f"   Center of mass: {ase_atoms.get_center_of_mass()}")
            print(f"   Total mass: {ase_atoms.get_masses().sum():.2f} amu")
            
            # Example ASE operations
            print("   ASE-specific operations:")
            print(f"   - Cell: {ase_atoms.get_cell()}")
            print(f"   - PBC: {ase_atoms.get_pbc()}")
            print(f"   - Symbols: {ase_atoms.get_chemical_symbols()}")
            
            # Calculate some properties
            if hasattr(ase_atoms, 'get_distances'):
                distances = ase_atoms.get_distances(0, range(len(ase_atoms)))
                print(f"   - Distances from first atom: {distances[:3]}...")
            
        except Exception as e:
            print(f"   ✗ Error with ASE: {e}")
    else:
        print("   Skipping ASE tests (not available)")
    print()
    
    # RDKit integration
    print("6. RDKit integration:")
    if rdkit_available:
        try:
            rdkit_mol = species.to_rdkit_mol()
            print(f"   ✓ Created RDKit Mol object")
            print(f"   Number of atoms: {rdkit_mol.GetNumAtoms()}")
            print(f"   Number of conformers: {rdkit_mol.GetNumConformers()}")
            
            # Example RDKit operations
            from rdkit import Chem
            try:
                smiles = Chem.MolToSmiles(rdkit_mol)
                print(f"   SMILES (explicit): {smiles}")
                
                # Convert to standard SMILES (implicit hydrogens)
                try:
                    standard_smiles = species.get_standard_smiles()
                    print(f"   SMILES (standard): {standard_smiles}")
                    if standard_smiles == 'CO' and len(species.get_xyz_data().atoms) == 6:
                        print("   ✓ Correctly identified as methanol (CH3OH)")
                except Exception as e:
                    print(f"   Standard SMILES conversion failed: {e}")
                
                if '.' in smiles:
                    print("   ⚠ Note: SMILES shows disconnected atoms - this may indicate")
                    print("     that bond inference failed for this molecular structure")
            except Exception as e:
                print(f"   SMILES generation failed: {e}")
            
            try:
                mw = Chem.rdMolDescriptors.CalcExactMolWt(rdkit_mol)
                print(f"   Molecular weight: {mw:.2f}")
            except Exception as e:
                print(f"   Molecular weight calculation failed: {e}")
            
            # Get conformer info
            try:
                conf = rdkit_mol.GetConformer(0)
                print(f"   Conformer 0 positions shape: {conf.GetPositions().shape}")
            except Exception as e:
                print(f"   Conformer access failed: {e}")
            
            # Calculate some descriptors
            print("   RDKit-specific operations:")
            try:
                print(f"   - Heavy atom count: {rdkit_mol.GetNumHeavyAtoms()}")
            except Exception as e:
                print(f"   - Heavy atom count failed: {e}")
            
            try:
                ring_count = Chem.rdMolDescriptors.CalcNumRings(rdkit_mol)
                print(f"   - Ring count: {ring_count}")
            except Exception as e:
                print(f"   - Ring count failed: {e}")
            
            try:
                aromatic_rings = Chem.rdMolDescriptors.CalcNumAromaticRings(rdkit_mol)
                print(f"   - Aromatic ring count: {aromatic_rings}")
            except Exception as e:
                print(f"   - Aromatic ring count failed: {e}")
            
        except Exception as e:
            print(f"   ✗ Error with RDKit: {e}")
            print("   Note: RDKit can be sensitive to molecular structure and may fail")
            print("   for some molecules. This is a known limitation of RDKit.")
            print("   Alternative: Use ASE or pymatgen for more robust molecular handling.")
    
    # Note about deprecation warnings
    if rdkit_available:
        print("   Note: You may see RDKit deprecation warnings about GetValence().")
        print("   These are harmless and come from RDKit's internal code.")
    else:
        print("   Skipping RDKit tests (not available)")
    print()
    
    # pymatgen integration
    print("7. pymatgen integration:")
    if pymatgen_available:
        try:
            pymatgen_mol = species.to_pymatgen_molecule()
            print(f"   ✓ Created pymatgen Molecule object")
            print(f"   Chemical formula: {pymatgen_mol.composition.formula}")
            print(f"   Number of sites: {len(pymatgen_mol.sites)}")
            print(f"   Center of mass: {pymatgen_mol.center_of_mass}")
            
            # Example pymatgen operations
            print("   pymatgen-specific operations:")
            print(f"   - Composition: {pymatgen_mol.composition}")
            print(f"   - Charge: {pymatgen_mol.charge}")
            print(f"   - Spin multiplicity: {pymatgen_mol.spin_multiplicity}")
            
            # Calculate some properties
            if hasattr(pymatgen_mol, 'get_distance'):
                if len(pymatgen_mol.sites) >= 2:
                    dist = pymatgen_mol.get_distance(0, 1)
                    print(f"   - Distance between first two atoms: {dist:.3f} Å")
            
        except Exception as e:
            print(f"   ✗ Error with pymatgen: {e}")
    else:
        print("   Skipping pymatgen tests (not available)")
    print()
    
    # File operations
    print("8. File operations:")
    try:
        species.save_xyz_to_file("species_from_api.xyz")
        print("   ✓ Saved XYZ data to species_from_api.xyz")
        
        # Verify file content
        with open("species_from_api.xyz", 'r') as f:
            file_content = f.read()
        print(f"   File size: {len(file_content)} characters")
        print(f"   First few lines:")
        lines = file_content.split('\n')[:5]
        for line in lines:
            print(f"     {line}")
        
    except Exception as e:
        print(f"   ✗ Error with file operations: {e}")
    print()
    
    # Demonstrate error handling
    print("9. Error handling examples:")
    
    # Try with species that has no XYZ data
    print("   a) Species without XYZ data:")
    try:
        species_no_xyz = await get_species("67-56-1*0")  # Without expand_xyz
        xyz_data = species_no_xyz.get_xyz_data()
        if xyz_data is None:
            print("      ✓ Correctly detected no XYZ data")
        else:
            print("      ✗ Unexpected: found XYZ data")
    except Exception as e:
        print(f"      ✗ Error: {e}")
    
    # Try to convert to ASE without XYZ data
    print("   b) Converting to ASE without XYZ data:")
    try:
        ase_atoms = species_no_xyz.to_ase_atoms()
        print("      ✗ Unexpected: conversion succeeded")
    except ValueError as e:
        print(f"      ✓ Correctly caught error: {e}")
    except Exception as e:
        print(f"      ✗ Unexpected error: {e}")
    print()
    
    # Clean up
    if os.path.exists("species_from_api.xyz"):
        os.remove("species_from_api.xyz")
        print("10. Cleaned up temporary files")
    
    print("\n=== Example completed ===")
    print("\nTo use these XYZ helpers in your own code:")
    print("1. Fetch species with: species = await get_species(atct_id, expand_xyz=True)")
    print("2. Convert to ASE: ase_atoms = species.to_ase_atoms()")
    print("3. Convert to RDKit: rdkit_mol = species.to_rdkit_mol()")
    print("4. Convert to pymatgen: pymatgen_mol = species.to_pymatgen_molecule()")
    print("5. Save to file: species.save_xyz_to_file('molecule.xyz')")

if __name__ == "__main__":
    asyncio.run(main())
