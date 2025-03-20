/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: Jan. 11, 2022                                                *
 * -------------------------------------------------------------------------- */

#pragma once
#include "Tools.h"
#include "Vec3.h"

/**
 * This class represents StructureVec3 of a system, including atom information,
 * positions, velocities, forces, and box vectors, which can be used to load
 * Structure file and store StructureVec3 data and save Structure/trajectory file.
 * The data are stored in OpenMM format.
 * Currently, only Gromacs Structure/trajectory files are supported.
 */
class StructureVec3 {
public:
    /**
    * Create a StructureVec3 object.
    *
    * This will initialize boxVectors and natoms to zero, and the other data, such
    * as positions, are empty. You can initialize the data using loadStructureVec3()
    * from a file, or using setNumAtoms(), setXXX(),..., to do this manully.
    */
    StructureVec3();
    ~StructureVec3() {}
    /**
    * Load data from Structure file.
    * This will override the existing StructureVec3 data (if any).
    *
    * Currently, the supported file is Gromacs gro file.
    *
    * @param file  the file name of structure file to load, the extension name
    *              of it indicates the format.
    */
    void loadStructure(const std::string& file);
    /**
    * Save data to StructureVec3 file.
    * This will override the existing file with same name (if any).
    *
    * Currently, the supported file is Gromacs gro file.
    *
    * @param file  the file name of StructureVec3 file to save, the extension name
    *              of it indicates the format.
    * @param step  the current MD simulation step
    * @param time  the current MD simulation in ps
    */
    void saveStructure(const std::string& file, int step, double time) const;
    /**
    * Save data to trajectory file.
    * When calling this function at first time, it will create a new file to
    * write trajectory (if an old file with the same name already exists, it
    * will be overwritten); When calling this function next time, the trajectory
    * will be appended to the file.
    *
    * If the velocities and/or forces are available and the format of trjectory
    * is supported to store, then, they will be saved to the file, too.
    *
    * Currently, the supported files are Gromacs gro/xtc/trr file under single
    * precision.
    *
    * @param file       the file name of trajectory file, the extension name of
    *                   it indicates the format.
    * @param step       the current MD simulation step
    * @param time       the current MD simulation in ps
    * @param group      the group to save to tajectory file. Currently, only atom
    *                   list can be accepted, such as "0-206,300". By default,
    *                   it is empty, means whole system will be written.
    * @param precision  the precision to save, typical and default value is 4,
    *                   which means 4 significant figures (0.001 nm). The allowed
    *                   value are 4, 5, 6. [xtc format only]
    */
    void saveTrajectory(const std::string& file, int step, double time, const std::string& group = "", int precision = 4) const;
    /**
     * Get the number of atoms of this structure.
     * If it is zero, which means the stucture data is not available.
     *
     * @return the number of atoms
     */
    int getNumAtoms() const {
        return natoms;
    }
    /**
     * Get reference to the atom information of each particle. [readonly]
     */
    const std::vector<std::string>& getAtomInfo() const {
        return atominfo;
    }
    /**
     * Get reference to the position (in nm) of each particle. [readonly]
     */
    const std::vector<Vec3>& getPositions() const {
        return positions;
    }
    /**
     * Get reference to the velocity (in nm/ps) of each particle. [readonly]
     */
    const std::vector<Vec3>& getVelocities() const {
        return velocities;
    }
    /**
     * Get reference to the force (in kj/mol/nm) of each particle. [readonly]
     */
    const std::vector<Vec3>& getForces() const {
        return forces;
    }
    /**
     * Get the vectors defining the axes of the periodic box (measured in nm).
     *
     * @param a      the vector defining the first edge of the periodic box
     * @param b      the vector defining the second edge of the periodic box
     * @param c      the vector defining the third edge of the periodic box
     */
    void getBoxVectors(Vec3& a, Vec3& b, Vec3& c) const;
    /**
     * Get reference to the atom information of each particle. [writable]
     */
    std::vector<std::string>& getAtomInfo() {
        return atominfo;
    }
    /**
     * Get reference to the position (in nm) of each particle. [writable]
     */
    std::vector<Vec3>& getPositions() {
        return positions;
    }
    /**
     * Get reference to the velocity (in nm/ps) of each particle. [writable]
     */
    std::vector<Vec3>& getVelocities() {
        return velocities;
    }
    /**
     * Get reference to the force (in kj/mol/nm) of each particle. [writable]
     */
    std::vector<Vec3>& getForces() {
        return forces;
    }
    /**
     * Get reference to the boxVectors (in nm). [writable]
     */
    auto getBoxVectors() -> Vec3(&)[3] {
        return boxVectors;
    }
    /**
     * Set the number of atoms of this structure to given number.
     *
     * It is only works for a empty structure whose natoms is zero.
     *
     * @param natoms the number of atoms
     */
    void setNumAtoms(int natoms);
    /**
     * Set the atom information of each particle.
     *
     * @param atominfo a vector whose length equals the number of particles
     */
    void setAtomInfo(const std::vector<std::string>& atominfo);
    /**
     * Set the positions of all particles (measured in nm).
     *
     * @param positions a vector whose length equals the number of particles
     */
    void setPositions(const std::vector<Vec3>& positions);
    /**
     * Set the velocities of all particles (measured in nm/ps).
     *
     * @param velocities a vector whose length equals the number of particles
     */
    void setVelocities(const std::vector<Vec3>& velocities);
    /**
     * Set the forces of all particles (measured in kj/mol/nm).
     *
     * @param forces a vector whose length equals the number of particles
     */
    void setForces(const std::vector<Vec3>& forces);
    /**
     * Set the vectors defining the axes of the periodic box (measured in nm).
     * The validity of vectors will be checked.
     *
     * @param a      the vector defining the first edge of the periodic box
     * @param b      the vector defining the second edge of the periodic box
     * @param c      the vector defining the third edge of the periodic box
     */
    void setBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c);

private:
    /**
     * Load atominfo, positions, velocities and box vectors from Gromacs structure file.
     *
     * Format Details of Gromacs Structure File (.gro):
     * 1st line: title string, optional time in ps after 't='
     * 2nd line: number of atoms (free format integer)
     * 3rd ~   : one line for each atom (fixed format): all columns are in a fixed
     * position, format will then be n+5 positions with n decimal places
     * (n+1 for velocities) in stead of 8 with 3 (with 4 for velocities), columns
     * contain the following information (from left to right):
     * residue index (5 positions, integer, from 1), residue name (5 characters),
     * atom name (5 characters), atom index (5 positions, integer),
     * position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
     * velocity (in nm/ps, x y z in 3 columns, each 8 positions with 4 decimal places)
     * Note that separate molecules or ions (e.g. water or Cl-) are regarded as residues,
     * and some fields may be written without spaces, (mind that when reading), for
     * example, there is always no space between residue index and residue name.
     * C format for writing: "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
     * After each atom line: box vectors (free format, space separated double),
     * values: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y),
     * the last 6 values may be omitted (they will be set to zero), the origin is
     * (0, 0, 0). Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0.
     *
     * Note:
     * 1. since there are position limitation for each filed, sometimes, there
     *    is no any space between two filed, for example, when atomindex is over
     *    99999, only last 5 digitals are shown and no space between the atomname
     *    and atomindex.
     * 2. If triclinic box is used, the box vectors loaded from Gromacs structure
     *    file be converted to OpenMM reduced form.
     * 3. If the file contains multi-frames, only the first frame will be loaded.
     *
     * @param file   the file name of Gromacs structure file
     */
    void readGROFile(const std::string& file);
    /**
     * Load atominfo, positions, velocities and box vectors from Gromacs structure file.
     *
     * Format Details of Gromacs Structure File (.ngro):
     * 1st line: title string, optional time in ps after 't='
     * 2nd line: number of atoms (free format integer)
     * 3rd ~   : one line for each atom (fixed format): all columns are in a fixed
     * position, format will then be n+5 positions with n decimal places
     * (n+1 for velocities) in stead of 8 with 3 (with 4 for velocities), columns
     * contain the following information (from left to right):
     * residue index (5 positions, integer, from 1), residue name (5 characters),
     * atom name (5 characters), atom index (5 positions, integer),
     * position (in nm, x y z in 3 columns, each 8 positions with 3 decimal places)
     * velocity (in nm/ps, x y z in 3 columns, each 8 positions with 4 decimal places)
     * Note that separate molecules or ions (e.g. water or Cl-) are regarded as residues,
     * and some fields may be written without spaces, (mind that when reading), for
     * example, there is always no space between residue index and residue name.
     * C format for writing: "%5d%-5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f"
     * After each atom line: box vectors (free format, space separated double),
     * values: v1(x) v2(y) v3(z) v1(y) v1(z) v2(x) v2(z) v3(x) v3(y),
     * the last 6 values may be omitted (they will be set to zero), the origin is
     * (0, 0, 0). Gromacs only supports boxes with v1(y)=v1(z)=v2(z)=0.
     *
     * Note:
     * 1. since there are position limitation for each filed, sometimes, there
     *    is no any space between two filed, for example, when atomindex is over
     *    99999, only last 5 digitals are shown and no space between the atomname
     *    and atomindex.
     * 2. If triclinic box is used, the box vectors loaded from Gromacs structure
     *    file be converted to OpenMM reduced form.
     * 3. If the file contains multi-frames, only the first frame will be loaded.
     *
     * @param file   the file name of Gromacs structure file
     */
    void readNGROFile(const std::string& file);
    /**
    * Save structure (single frame) or trajectory (multi frames) to Gromacs gro
    * file (ASCII, positions/velocities, fixed precision = 4 significant figures).
    *
    * If velocities are available, then they will be saved, too.
    *
    * @param file    the file name of gro file
    * @param step    the current MD simulation step
    * @param time    the current MD simulation in ps
    * @param indices the atom indcies to save to tajectory file. By default,
    *                it is empty, means whole system will be written.
    * @param isTraj  By default, it is false, which means to save structure
    *                file (single frame). If it is true, then next frame will be
    *                appended to the file, as trajectory file.
    */
    void writeGROFile(const std::string& file, int step, double time, const std::vector<int>& indices = std::vector<int>(), bool isTraj = false) const;
    /**
    * Save trajectory to Gromacs xtc file (compressed binary, positions only,
    * any precision) using Gromacs XDR library.
    *
    * The XDR library can be downloaded from the Gromacs FTP site:
    * ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.4.tar.gz
    * The installation and usage of it, please refer to the source file directly.
    * Note: the XDR library used here has been adsorbed in the chemfiles library,
    * and the header files in the XDR library are included in "chemfiles.hpp".
    * See the introduction of chemfiles in https://github.com/chemfiles/chemfiles
    *
    * @param file      the file name of xtc file
    * @param step      the current MD simulation step
    * @param time      the current MD simulation in ps
    * @param indices   the atom indcies to save to tajectory file. By default,
    *                  it is empty, means whole system will be written.
    * @param precision the precision to save, typical and default value is 4,
    *                  which means 4 significant figures (0.001 nm). The allowed
    *                  value are 4, 5, 6.
    */
    void writeXTCFile(const std::string& file, int step, double time, const std::vector<int>& indices = std::vector<int>(), int precision = 4) const;
    /**
    * Read trajectory from Gromacs trr file (binary, positions/velocities/forces,
    * full precision) using Gromacs XDR library.
    *
    * The XDR library can be downloaded from the Gromacs FTP site:
    * ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.4.tar.gz
    * The installation and usage of it, please refer to the source file directly.
    * Note: the XDR library used here has been adsorbed in the chemfiles library,
    * and the header files in the XDR library are included in "chemfiles.hpp".
    * See the introduction of chemfiles in https://github.com/chemfiles/chemfiles
    *
    * If velocities and/or forces are available, then they will be loaded, too.
    *
    * @param file    the file name of trr file
    * @param step    the MD simulation step to load
    * @param time    the MD simulation time in ps to load
    */
    void readTRRFile(const std::string& file, int step, double time);
    /**
    * Save trajectory to Gromacs trr file (binary, positions/velocities/forces,
    * full precision) using Gromacs XDR library.
    *
    * The XDR library can be downloaded from the Gromacs FTP site:
    * ftp://ftp.gromacs.org/pub/contrib/xdrfile-1.1.4.tar.gz
    * The installation and usage of it, please refer to the source file directly.
    * Note: the XDR library used here has been adsorbed in the chemfiles library,
    * and the header files in the XDR library are included in "chemfiles.hpp".
    * See the introduction of chemfiles in https://github.com/chemfiles/chemfiles
    *
    * If velocities and/or forces are available, then they will be saved, too.
    *
    * @param file    the file name of trr file
    * @param step    the current MD simulation step
    * @param time    the current MD simulation in ps
    * @param indices the atom indcies to save to tajectory file. By default,
    *                it is empty, means whole system will be written.
    */
    void writeTRRFile(const std::string& file, int step, double time, const std::vector<int>& indices = std::vector<int>()) const;
    /**
     * Check the vectors defining the axes of the periodic box (measured in nm).
     *
     * Triclinic boxes are supported, but the vectors must satisfy certain requirements.
     * In particular, a must point in the x direction, b must point "mostly" in the
     * y direction, and c must point "mostly" in the z direction.
     * See the OpenMM documentation for details.
     *
     * If the box vectors are not satisfied the requirements, an exception will
     * be thrown. If they are OK, then do nothing.
     *
     * @param a    the vector defining the first edge of the periodic box
     * @param b    the vector defining the second edge of the periodic box
     * @param c    the vector defining the third edge of the periodic box
     */
    void checkBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c) const;

private:
    // The numbers of atoms or particles for this structure. When creating a
    // object, the initial value is zero, and it will be setted automatically
    // after loading data from file or setting some data according to the length
    // of data. Or you can use setNumAtoms() to set it manully for an empty object.
    // Once it has been setted to non zero, it cannot be changed again,
    // all the subsequent changes of data should be consistent with it.
    int natoms;
    // atominfo represents the information of each atom from Gromacs gro file,
    // includes resid, resname, atomname and atomindex, which can be loaded from
    // Gromacs gro file, and will be used in the saving of Gromacs gro file.
    std::vector<std::string> atominfo;
    // The coordinate (in nm) of each atom in Vec3 format.
    std::vector<Vec3> positions;
    // The velocity in (nm/ps) of each atom in Vec3 format.
    std::vector<Vec3> velocities;
    // The forces in (kj/mol/nm) of each atom in Vec3 format.
    std::vector<Vec3> forces;
    // The three box vectors in OpenMM reduced form.
    // For rectangular/cubic box, it has same form as Gromacs, and no conversion
    // is required. While for triclinic box, although the requiremnt is same as
    // for Gromacs and OpenMM form, but some transformation may be needed to make
    // sure they're in the reduced form required by OpenMM when apply it to
    // OpenMM System object. However, using OpenMM form directly is good when
    // output box vectors to structrue or trajectory file since it can be
    // recognized by VMD and Gromacs [has been tested by zhubin].
    // If periodic boundary conditions (PBC) is not used for this system, all
    // the box vectors are 0.
    Vec3 boxVectors[3];
};