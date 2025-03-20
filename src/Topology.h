/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Feb. 13, 2021                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "Tools.h"

/**
 * This class represents the force filed parameters of one system loaded from
 * topology file, such as Gromacs topology file (.top). The definition of a
 * Topology involves four members:
 * 1. The ForceFiled are used to describe the type of force filed used. But, only
 *    AMBER/GAFF force filed is supported now.
 * 2. The AtomTypes stores the name, mass, and van der Waals parameters of atom.
 * 3. The Melecules records the name, type, and number of molecule of this system.
 * 4. The MoleculeTypes stores the atoms (charge), and bonded parameters (such as
 *    bond, angle, dihedral, and so on) of one kind of molecule.
 *
 * Note:
 * 1. this class is designed according to the format of Gromacs topology
 *    file, and the units are that used in Gromacs, i.e., distance in nm, energy in
 *    kj/mol, angle/dihedral in degree, force constant in kj/mol/nm^2 or kj/mol/rad^2.
 *    The unit for angle/dihedral in OpenMM is radian, therefore a conversion from
 *    degree to radian for angle is requested when creating an OpenMM System object.
 * 2. Only only a standalone Gromacs topology file is supported now, which means
 *    the description of "#include" or "#define" can not be processed now.
 *    To create this standalone topology file, ParmEd can be used to convert the
 *    Amber .prmtop and .inpcrd or .rst files into Gromacs .top and .gro files.
 */
class Topology {
    /**
     * Internal class of Topology used to describe the type of force filed used. The
     * definition of ForceFiled involves the following members:
     * 1. The name of force filed.
     * 2. The form of the repulsion and dispersion term of van der Walls interaction.
     *    It can be either the (1) Lennard-Jones (or 6-12 interaction) or the (2)
     *    Buckingham (or exp-6 potential) potential. Buckingham is more accurate but
     *    expensive and seldom used.
     *    Only Lennard-Jones potential is supported now.
     * 3. The combination rule to construct the parameter matrix for the nonbonded
     *    Lennard-Jones parameters. It can be (1) used for GROMOS force filed. (2)
     *    Lorentz-Berthelot rule, arithmetic average for sigma, and geometric average
     *    for epsilon, i.e., sigma=(sigma1+siam2)/2, epsilon=sqrt(epsilon1+epsilon2)
     *    which is used for AMBER or CHARMM force filed. (3) geometric average for
     *    both parameter, used for OPLS force filed.
     *    Only Lorentz-Berthelot rule is supported now.
     * 4. Whether or not generating vdW parameters for 1-4 interaction automatically.
     *    For GROMOS force field, should be false, since these paramters have been
     *    defined in its topology file ([ pairtypes ]). For AMBER/CHARMM/OPLS-AA force
     *    field, it should be true.
     *    This parameter is not used now.
     * 5. The Scale factors used for 1-4 Coulomb and Lennard Jones interactions.
     *    For AMBER force field, they should be 0.833333 and 0.5, respectively.
     *    For GROMOS or CHARMM force filed, both of them should be 1.0.
     *    For OPLS-AA force field, both of them should be 0.5.
     * Note that, only AMBER/GAFF force filed is supported now.
     * @private
     */
    struct ForceFiled { std::string name; int vdWType, combinationRule; bool genPairs; double Coulomb14Scale, LennardJones14Scale; };
    /**
     * Internal class of Topology used to records the name, mass, and van der Waals
     * parameters (sigma/nm and epsilon/kJ·mol-1) of atom in a system. sigma is the
     * distance at which the energy equals zero and epsilon sets the strength of the
     * interaction.
     * @private
     */
    struct AtomTypes { std::string atomType; int atomicNumber; double mass, sigma, epsilon; };
    /**
     * Internal class of Topology used to records the name, type and number of molecule
     * in a system.
     * Note that, the order of molecules listed here must be same with that in structure file.
     * And the name of molecules must be same with that in MoleculeTypes. But, one
     * molecule type can be defined several times in molecules to follow the order
     * of molecules in structure file. Don't include space in the molecule name.
     * @private
     */
    struct Molecules { std::string moleculeName; int moleculeNumber, moleculeTypeIndex; };
    /**
     * Internal class of Topology used to create one entry of molecule type, which
     * stores the atoms (charge), and bonded parameters (such as bond, angle, dihedral
     * , and so on) of one kind of molecule.
     * @private
     */
public:
    struct MoleculeTypes {
        /**
         * Internal class of MoleculeTypes used to records the parameters of the type,
         * and charge of atoms in one molecule.
         * Here, the atomTypeIndex records the index of atomType in AtomTypes.
         * @private
         */
        struct Atoms { int atomTypeIndex; std::string atomType; double charge; };
        /**
         * Internal class of MoleculeTypes used to records the parameters of the Harmonic
         * bond stretching (2-body) interaction. The form for Harmonic bond potential
         * is: E = 1/2*k(x-x0)^2, in which k is force constant (in kj/mol/nm^2), and x0
         * is equilibrium distance (in nm). Note that the atom index counts from 1.
         * Here, HydrogenTag records which atom is Hydrogen: 0: no-H; 1: atomI; 2: atomJ.
         * @private
         */
        struct Bonds { int atomIndex[2], HydrogenTag; double distance, forceConstant; };
        /**
         * Internal class of MoleculeTypes used to records the parameters of 1-4 (2-body)
         * interaction which should be calculated by multiplying the scale factor for Coulomb
         * and Lennard Jones interactions, for example, the scale factors used for Coulomb
         * and Lennard Jones interactions by AMBER/GAFF force filed are 0.5 and 0.833333,
         * respectively. And they should be excluded from the normal nonbonded interactions.
         * @private
         */
        struct Pairs { int atomIndex[2]; };
        /**
         * Internal class of MoleculeTypes used to records the parameters of the Harmonic
         * bond angle (3-body) interaction. The form for Harmonic angle potential is:
         * E = 1/2*k(c-c0)^2, in which k is force constant (in kj/mol/rad^2), and c0 is
         * equilibrium angle (in degree). Here, isHAngle means H-X-H or H-O-X (where X
         * is an arbitrary atom), which can be constrainted.
         * @private
         */
        struct Angles { int atomIndex[3]; double angle, forceConstant; bool isHAngle; };
        /**
         * Internal class of MoleculeTypes used to records the rigid water parameters,
         * includes the distance (in nm) between O-H and H-H distances. SETTLE algorithm
         * is an analytical solution of SHAKE, specifically for water.
         * Note that the parameters of water model converted from AMBER can only be used
         * as rigid water. Since the force constant for flexible water is not right.
         * But, you can copy these parameters from Gromacs file, such as tip3p.itp,
         * spc.itp, if flexible water model is used.
         * @private
         */
        // TODO: support 4-point water model.
        struct Settles { double distanceOH, distanceHH; };
        /**
         * Internal class of MoleculeTypes used to records the parameters of periodic type
         * dihedral (4-body) interaction. The form for periodic type dihedral potential is:
         * E = k(1+cos(n*c-c0)), where c is the dihedral angle (in degree) formed by the
         * four particles, c0 is the phase offset, n is the periodicity, and k is the
         * force constant (in kj/mol/rad^2).
         * Note that in AMBER/GAFF proper and improper dihedral angles are defined by
         * the same functional form (periodic type).
         * @private
         */
        struct Dihedrals { int atomIndex[4], periodicity; double dihedral, forceConstant; };

        /**
         * Internal class of MoleculeTypes used to records the parameters of extra
         * exclusions within a molecule. It starts with one atom index, followed by one
         * or more atom indices. All non-bonded interactions between the first atom and
         * the other atoms will be excluded for generating exclusions.
         * @private
         */
        struct Exclusions { std::vector<int> atomIndex; };
        /**
         * Internal class of MoleculeTypes used to records the parameters of extra
         * constraints for bonds with or without chemical bond, the former will be used
         * for generating exclusions.
         * @private
         */
        struct Constraints { int atomIndex[2]; bool useGenExclusion; double distance; };
        /**
         * Internal class of MoleculeTypes used to records the parameters of the type,
         * and polarizability of atoms in one molecule.
         * Here, the atomTypeIndex records the index of atomType in AtomTypes.
         * @private
         */
        struct Polars { int atomTypeIndex; std::string atomType; double polarizability; };
        /**
         * Internal class of MoleculeTypes used to records the parameters of the type,
         * and new vdw model parameters of atoms in one molecule.
         * Here, the atomTypeIndex records the index of atomType in AtomTypes.
         * @private
         */
        struct NewVDW { int atomTypeIndex; double vdwparam0; double vdwparam1; 
                        double vdwparam2; double vdwparam3; double vdwparam4; double vdwparam5; double vdwparam6;
                        double vdwparam7; double vdwparam8;};
        
        std::string              moleculeName;
        std::vector<Atoms>       atoms;
        std::vector<Bonds>       bonds;
        std::vector<Pairs>       pairs;
        std::vector<Angles>      angles;
        std::vector<Settles>     settles;
        std::vector<Dihedrals>   dihedrals;
        std::vector<Exclusions>  exclusions;
        std::vector<Constraints> constraints;
        std::vector<Polars>      polars;
        std::vector<NewVDW>      newvdw;
    };
public:
    std::string                systemName;
    ForceFiled                 forceFiled;
    std::vector<AtomTypes>     atomTypes;
    std::vector<Molecules>     molecules;
    std::vector<MoleculeTypes> moleculeTypes;
public:
    /**
    * Create a Topology object.
    *
    * After creating a Topology object, you can call loadStructureFile() to
    * load data from a file.
    */
    Topology() {}
    ~Topology() {}
    /**
    * Load data from topology file.
    *
    * The supported file is Gromacs standalone topology file.
    *
    * @param file  the file name of topology file
    */
    void loadTopologyFile(const std::string& file);
    /**
     * Compute and return total physical mass (in u) of system based on original
     * topology, which can be used to compute the density.
     *
     * Unit of atomic mass is amu/u or Dalton, which is 1/12 of C-12 mass.
     * 1 u = 1.660538921e-27 kg. (Reference: Gromacs manual)
     *
     * @return    the total physical mass (in u)
     */
    int computeSystemMass() const;
    /**
    * Print the information of topology, such as force filed type, the molecules
    * in this system, and so on.
    */
    void printTopologyInfo() const;
private:
    /**
     * Load force filed paramters from Gromacs topology file.
     *
     * @param topfile   the file stream of Gromacs topology file
     */
    void loadGromacsTopology(std::ifstream& topfile);
};