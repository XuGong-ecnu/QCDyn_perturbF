/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Xiaofang Zhang @Sun Group @NYU-SH                                       *
 * Last updated: Dec. 15, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "StructureVec3.h"
#include "chemfiles.hpp"

StructureVec3::StructureVec3() {
    natoms = 0;
    for (int i = 0; i < 3; ++i)
        boxVectors[i] = Vec3(0, 0, 0);
}

void StructureVec3::loadStructure(const std::string& file) {
    const std::string suffix = GetFileSuffix(file); // extension name
    if (suffix == "gro")
        {
           readGROFile(file);
        }
    else if (suffix == "inpcrd")
        throw std::runtime_error("ERROR: The Amber inpcrd structure file is not yet supported.");
    else if (suffix == "pdb")
        throw std::runtime_error("ERROR: The PDB structure file is not yet supported.");
    // add the suffix == "ngro"
    else if (suffix == "ngro")
        readNGROFile(file);
    else
        throw std::runtime_error("ERROR: Unkown structure file format: " + suffix);
}

void StructureVec3::saveStructure(const std::string& file, int step, double time) const {
    if (natoms == 0)
        throw std::runtime_error("ERROR: Called saveStructure() on a Structure with empty positions.");
    if (this->positions.size() != natoms)
        throw std::runtime_error("ERROR: Called saveStructure() on a Structure with wrong number of positions.");
    const std::string suffix = GetFileSuffix(file); // extension name
    if (suffix == "gro")
        writeGROFile(file, step, time);
    else if (suffix == "pdb")
        throw std::runtime_error("ERROR: The PDB structure file is not yet supported.");
    else
        throw std::runtime_error("ERROR: Unkown structure file format: " + suffix);
}

void StructureVec3::saveTrajectory(const std::string& file, int step, double time, const std::string& group, int precision) const {
    if (natoms == 0)
        throw std::runtime_error("ERROR: Called saveTrajectory() on a Structure with empty positions.");
    if (this->positions.size() != natoms)
        throw std::runtime_error("ERROR: Called saveTrajectory() on a Structure with wrong number of positions.");
    // Decide which atom groups shold be wriiten to trajectory fille.
    static int count = 0;
    static std::vector<int> indices;
    // TODO: more advance selection
    if (count == 0) { // do this one time for one system at the first time
        // if  group is empty or system, means whole system to save.
        if (group != "" && group != "system") {
            GetIndexList(indices, group); // Get indices according group selections.
            if (indices.size() > natoms)
                throw std::runtime_error("ERROR: The atom indices in the group to save "
                    "trajectory file is out of range.");
        }
        count++;
    }
    const std::string suffix = GetFileSuffix(file); // extension name
    if (suffix == "gro")
        writeGROFile(file, step, time, indices, true);
    else if (suffix == "xtc")
        writeXTCFile(file, step, time, indices, precision);
    else if (suffix == "trr")
        writeTRRFile(file, step, time, indices);
    else if (suffix == "tng")
        throw std::runtime_error("ERROR: The Gromacs TNG trajectory file is not yet supported.");
    else if (suffix == "crd" || suffix == "mdcrd")
        throw std::runtime_error("ERROR: The Amber ASCII trajectory file is not yet supported.");
    else if (suffix == "nc")
        throw std::runtime_error("ERROR: The Amber NetCDF trajectory file is not yet supported.");
    else if (suffix == "pdb")
        throw std::runtime_error("ERROR: The PDB trajectory file is not yet supported.");
    else
        throw std::runtime_error("ERROR: Unkown trajectory file format: " + suffix);
}

void StructureVec3::getBoxVectors(Vec3& a, Vec3& b, Vec3& c) const {
    a = this->boxVectors[0];
    b = this->boxVectors[1];
    c = this->boxVectors[2];
}

void StructureVec3::setNumAtoms(int natoms) {
    if (this->natoms != 0)
        throw std::runtime_error("ERROR: You cannot change the number of atoms of a existing structure.");
    if (natoms < 0)
        throw std::runtime_error("ERROR: Called setNumAtoms() with a negative number of atoms.");
    this->natoms = natoms;
}

void StructureVec3::setAtomInfo(const std::vector<std::string>& atominfo) {
    if (natoms != 0 && atominfo.size() != natoms)
        throw std::runtime_error("ERROR: Called setAtomInfo() on a Structure with the wrong number of atominfo.");
    else if (natoms == 0)
        natoms = positions.size();
    this->atominfo = atominfo;
}

void StructureVec3::setPositions(const std::vector<Vec3>& positions) {
    if (natoms != 0 && positions.size() != natoms)
        throw std::runtime_error("ERROR: Called setPositions() on a Structure with the wrong number of positions.");
    else if (natoms == 0)
        natoms = positions.size();
    this->positions = positions;
}

void StructureVec3::setVelocities(const std::vector<Vec3>& velocities) {
    if (natoms != 0 && velocities.size() != natoms)
        throw std::runtime_error("ERROR: Called setVelocities() on a Structure with the wrong number of velocities.");
    else if (natoms == 0)
        natoms = velocities.size();
    this->velocities = velocities;
}

void StructureVec3::setForces(const std::vector<Vec3>& forces) {
    if (natoms != 0 && forces.size() != natoms)
        throw std::runtime_error("ERROR: Called setForces() on a Structure with the wrong number of forces.");
    else if (natoms == 0)
        natoms = forces.size();
    this->forces = forces;
}

void StructureVec3::setBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c) {
    checkBoxVectors(a, b, c);
    this->boxVectors[0] = a;
    this->boxVectors[1] = b;
    this->boxVectors[2] = c;
}

void StructureVec3::checkBoxVectors(const Vec3& a, const Vec3& b, const Vec3& c) const {
    if (a[1] != 0.0 || a[2] != 0.0)
        throw std::runtime_error("ERROR: First periodic box vector must be parallel to x.");
    if (b[2] != 0.0)
        throw std::runtime_error("ERROR: Second periodic box vector must be in the x-y plane.");
    if (a[0] <= 0.0 || b[1] <= 0.0 || c[2] <= 0.0 || a[0] < 2*fabs(b[0]) || a[0] < 2*fabs(c[0]) || b[1] < 2*fabs(c[1]))
        throw std::runtime_error("ERROR: Periodic box vectors must be in reduced form.");
}

void StructureVec3::readGROFile(const std::string& file) {
    if (GetFileSuffix(file) != "gro")
        throw std::runtime_error("ERROR: Called readGROFile() with wrong extension name.");
    std::ifstream grofile(file.c_str());
    if (!grofile)
        throw std::runtime_error("ERROR: Can't open structure file: " + file);
    int lineNumber = 0;
    bool hasVel = false; // true if this file contains velocities
    for (std::string line; getline(grofile, line); ) {
        lineNumber++; // Let line number start from 1.
        // First line is comment.
        if (lineNumber == 1) continue;
        // Get the total number of atoms at line 2.
        else if (lineNumber == 2) {
            int natoms_gro = std::stoi(line);
            if (natoms != 0 && natoms_gro != natoms)
                throw std::runtime_error("ERROR: Called readGROFile() but got inconsistent natoms from file: " + file);
            else if (natoms == 0)
                natoms = natoms_gro;
            atominfo.resize(natoms);
            positions.resize(natoms);
            // Don't know if this file contains velocities, clear it firstly
            velocities.clear();
            forces.clear(); // No forces data in gro file.
        }
        // Line 3 - Line natoms+2: atominfo, positions and/or velocities
        else if (lineNumber > 2 && lineNumber < (natoms + 3)) {
            int numChars = line.size();
            if (numChars != 44 && numChars != 68)
                throw std::runtime_error("ERROR: Called readGROFile() but found illegal "
                    "data at line " + std::to_string(lineNumber) + " of file: " + file);
            if (lineNumber == 3 && numChars == 68) { // Found data of velocities
                hasVel = true;                       // do this one time only
                velocities.resize(natoms);
            }
            std::string atom = line.substr(0, 20);
            double px = std::stod(line.substr(20, 8));
            double py = std::stod(line.substr(28, 8));
            double pz = std::stod(line.substr(36, 8));
            atominfo[lineNumber-3] = atom;
            positions[lineNumber-3] = Vec3(px, py, pz);
            if (hasVel) {
                double vx = std::stod(line.substr(44, 8));
                double vy = std::stod(line.substr(52, 8));
                double vz = std::stod(line.substr(60, 8));
                velocities[lineNumber-3] = Vec3(vx, vy, vz);
            }
        }
        // Line natoms+3 stores the box vectors. It could be a cubic/
        // rectangular box with 3 values, or a triclinic box with 9 values.
        else if (lineNumber == (natoms + 3)) {
            std::vector<double> values;
            SplitLine(values, line);
            Vec3 a, b, c;
            // This is cubic/rectangular box
            if (values.size() == 3) {
                a = Vec3(values[0], 0.0, 0.0);
                b = Vec3(0.0, values[1], 0.0);
                c = Vec3(0.0, 0.0, values[2]);
            }
            // Perhaps it is a triclinic box if the last 6 values are not zero.
            else if (values.size() == 9) {
                a = Vec3(values[0], values[3], values[4]);
                b = Vec3(values[5], values[1], values[6]);
                c = Vec3(values[7], values[8], values[2]);
                // Make sure they're in the reduced form required by OpenMM.
                // Reference: OpenMM python wrapper: unitcell.py
                c = c - b*round(c[1]/b[1]);
                c = c - a*round(c[0]/a[0]);
                b = b - a*round(b[0]/a[0]);
            }
            else
                throw std::runtime_error("ERROR: Called readGROFile() but found illegal value "
                    "for box vectors at line " + std::to_string(lineNumber) + " of file: " + file);
            // Check if the box vectors are valid. (to satisfy the OpenMM requirements)
            checkBoxVectors(a, b, c);
            this->boxVectors[0] = a;
            this->boxVectors[1] = b;
            this->boxVectors[2] = c;
            break; // if multi-frames in .gro file, only read 1st frame.
        }
    }
    grofile.close();
    // Check the number of data from this file.
    if (lineNumber != (natoms + 3))
        throw std::runtime_error("ERROR: Called readGROFile() but the number of "
            "positions or velocities doesn't math the number of atoms.");
}

void StructureVec3::readNGROFile(const std::string& file) {
    if (GetFileSuffix(file) != "ngro")
        throw std::runtime_error("ERROR: Called readGROFile() with wrong extension name.");
    std::ifstream ngrofile(file.c_str());
    if (!ngrofile)
        throw std::runtime_error("ERROR: Can't open structure file: " + file);
    int lineNumber = 0;
    bool hasVel = false; // true if this file contains velocities
    for (std::string line; getline(ngrofile, line); ) {
        lineNumber++; // Let line number start from 1.
        // First line is comment.
        if (lineNumber == 1) continue;
        // Get the total number of atoms at line 2.
        else if (lineNumber == 2) {
            int natoms_gro = std::stoi(line);
            if (natoms != 0 && natoms_gro != natoms)
                throw std::runtime_error("ERROR: Called readNGROFile() but got inconsistent natoms from file: " + file);
            else if (natoms == 0)
                natoms = natoms_gro;
            atominfo.resize(natoms);
            positions.resize(natoms);
            // Don't know if this file contains velocities, clear it firstly
            velocities.clear();
            forces.clear(); // No forces data in gro file.
        }
        // Line 3 - Line natoms+2: atominfo, positions and/or velocities
        else if (lineNumber > 2 && lineNumber < (natoms + 3)) {
            std::vector<std::string> fields;
            SplitLine(fields, line);
            if (fields.size() == 9)  { // Found data of velocities
                hasVel = true;                       // do this one time only
                velocities.resize(natoms);
            }
            std::string atom = line.substr(0, 20);
            double px = std::stod(fields[3]);
            double py = std::stod(fields[4]);
            double pz = std::stod(fields[5]);
            atominfo[lineNumber-3] = atom;
            positions[lineNumber-3] = Vec3(px, py, pz);
            if (hasVel) {
                double vx = std::stod(fields[6]);
                double vy = std::stod(fields[7]);
                double vz = std::stod(fields[8]);
                velocities[lineNumber-3] = Vec3(vx, vy, vz);
            }
        }
        // Line natoms+3 stores the box vectors. It could be a cubic/
        // rectangular box with 3 values, or a triclinic box with 9 values.
        else if (lineNumber == (natoms + 3)) {
            std::vector<double> values;
            SplitLine(values, line);
            Vec3 a, b, c;
            // This is cubic/rectangular box
            if (values.size() == 3) {
                a = Vec3(values[0], 0.0, 0.0);
                b = Vec3(0.0, values[1], 0.0);
                c = Vec3(0.0, 0.0, values[2]);
            }
            // Perhaps it is a triclinic box if the last 6 values are not zero.
            else if (values.size() == 9) {
                a = Vec3(values[0], values[3], values[4]);
                b = Vec3(values[5], values[1], values[6]);
                c = Vec3(values[7], values[8], values[2]);
                // Make sure they're in the reduced form required by OpenMM.
                // Reference: OpenMM python wrapper: unitcell.py
                c = c - b*round(c[1]/b[1]);
                c = c - a*round(c[0]/a[0]);
                b = b - a*round(b[0]/a[0]);
            }
            else
                throw std::runtime_error("ERROR: Called readNGROFile() but found illegal value "
                    "for box vectors at line " + std::to_string(lineNumber) + " of file: " + file);
            // Check if the box vectors are valid. (to satisfy the OpenMM requirements)
            checkBoxVectors(a, b, c);
            this->boxVectors[0] = a;
            this->boxVectors[1] = b;
            this->boxVectors[2] = c;
            break; // if multi-frames in .gro file, only read 1st frame.
        }
    }
    ngrofile.close();
    // Check the number of data from this file.
    if (lineNumber != (natoms + 3))
        throw std::runtime_error("ERROR: Called readNGROFile() but the number of "
            "positions or velocities doesn't math the number of atoms.");
}

void StructureVec3::writeGROFile(const std::string& file, int step, double time, const std::vector<int>& indices, bool isTraj) const {
    if (GetFileSuffix(file) != "gro")
        throw std::runtime_error("ERROR: Called writeGROFile() with wrong extension name.");
    static int frame = 0; // current frame index in trajectory file
    static std::string file_prev = file; // the filename of previous call
    if (file_prev != file) { // if it is not a same file,
        frame = 0;           // we reset the frame to 0,
        file_prev = file;    // and set current file name as previous one.
    }
    // the number of atoms to save.
    const int nindices = indices.size();
    if (nindices > this->natoms)
        throw std::runtime_error("ERROR: The atom indices in the group to save "
            "trajectory file is out of range.");
    const int natoms = nindices == 0 ? this->natoms : nindices;
    // write for frame = 0, and append for frame > 0
    FILE* groFile = CheckFile(file, isTraj ? frame : 0);
    // Line 1 is comment including current time (in ps).
    fprintf(groFile, "GROningen MAchine for Chemical Simulation t=%10.5f step= %d\n", time, step);
    // Line 2 is total atoms of system.
    fprintf(groFile, "%d\n", natoms);
    // Line from 3 is atomifno (resid, resname, atomname, atomindex), position, (and velocity).
    if (natoms == this->natoms) // whole system
        if (velocities.size() == natoms) // save velocities if available.
            for (int i = 0; i < natoms; ++i)
                fprintf(groFile, "%s%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", atominfo[i].c_str(),
                positions[i][0], positions[i][1], positions[i][2], velocities[i][0],
                velocities[i][1], velocities[i][2]);
        else // positions only
            for (int i = 0; i < natoms; ++i)
                fprintf(groFile, "%s%8.3f%8.3f%8.3f\n", atominfo[i].c_str(),
                    positions[i][0], positions[i][1], positions[i][2]);
    else // the selected atoms
        if (velocities.size() == natoms) // save velocities if available.
            for (int i = 0; i < natoms; ++i)
                fprintf(groFile, "%s%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f\n", atominfo[indices[i]].c_str(),
                positions[indices[i]][0], positions[indices[i]][1], positions[indices[i]][2],
                velocities[indices[i]][0], velocities[indices[i]][1], velocities[indices[i]][2]);
        else // positions only
            for (int i = 0; i < natoms; ++i)
                fprintf(groFile, "%s%8.3f%8.3f%8.3f\n", atominfo[i].c_str(),
                    positions[indices[i]][0], positions[indices[i]][1], positions[indices[i]][2]);
    // Next line is the box vectors.
    // The validity of box vectors will not be checked, since non-PBC simulation is
    // possible, in the case of non-PBC, the box vectors are zero.
    // If it is a rectangular/cubic box, only 3 diagonal values are needed.
    fprintf(groFile, "%10.5f%10.5f%10.5f", boxVectors[0][0], boxVectors[1][1], boxVectors[2][2]);
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j) {
            if (i == j) continue;
            if (boxVectors[i][j] != 0.0) { // triclinic box, has none-zero off-diagonal value
                fprintf(groFile, "%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f", boxVectors[0][1],
                    boxVectors[0][2],boxVectors[1][0], boxVectors[1][2], boxVectors[2][0], boxVectors[2][1]);
                goto ENDLOOP;
            }
        }
    ENDLOOP:
    fprintf(groFile, "\n");
    fclose(groFile);
    if (isTraj) // for multi-frame trajectory file
        frame++;
}

void StructureVec3::writeXTCFile(const std::string& file, int step, double time, const std::vector<int>& indices, int precision) const {
    if (GetFileSuffix(file) != "xtc")
        throw std::runtime_error("ERROR: Called writeXTCFile() with wrong extension name.");
    if (precision < 4 || precision > 6)
        throw std::runtime_error("ERROR: The precision of xtc trajectory file should be between 4 and 6.");
    static int frame = 0; // current frame index in trajectory file
    static std::string file_prev = file; // the filename of previous call
    if (file_prev != file) { // if it is not a same file,
        frame = 0;           // we reset the frame to 0,
        file_prev = file;    // and set current file name as previous one.
    }
    // the number of atoms to save, never change in one trajectory file
    const int nindices = indices.size();
    if (nindices > this->natoms)
        throw std::runtime_error("ERROR: The atom indices in the group to save "
            "trajectory file is out of range.");
    static const int natoms = nindices == 0 ? this->natoms : nindices;
    const std::string mode = frame == 0 ? "w" : "a"; // "w" for writing, "a" for append.
    XDRFILE* xd = xdrfile_open(file.c_str(), mode.c_str()); // Open a file for xdr I/O
    if (xd == NULL)
        throw std::runtime_error("ERROR: Error when opening file to write xtc file.");
    // Note: for triclinic box, the boxVectors here in OpenMM reduced form, but
    // no conversion is needed since it can be recognized by VMD and Gromacs.
    static float box[3][3]; // box matrix for xtc file
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            box[i][j] = this->boxVectors[i][j]; // for rectangular/cubix only
    // The coordinates themselves stored in reduced precision for xtc file.
    static std::vector<rvec> x(natoms); // rvec is float rvec[3]
    if (natoms == this->natoms) // whole system
        for (int i = 0; i < natoms; ++i)
            for (int j = 0; j < 3; ++j)
                x[i][j] = this->positions[i][j];
    else                    // the selected atoms
        for (int i = 0; i < natoms; ++i)
            for (int j = 0; j < 3; ++j)
                x[i][j] = this->positions[indices[i]][j];
    // Write or append a frame to xtc file. Return 0 if success.
    // Here, x.data() return a pointer to the first element, i.e. &x[0].
    // The precision for xtc, default is 1000, which means 4 significant figures.
    float prec = pow(10, precision);
    int isOK = write_xtc(xd, natoms, step, time, box, x.data(), prec);
    if (isOK != 0)
        throw std::runtime_error("ERROR: Error when writting data to xtc file.");
    xdrfile_close(xd); // Close the file for xdr I/O
    frame++;
}

// TODO: not finished
void StructureVec3::readTRRFile(const std::string& file, int step, double time) {
    if (GetFileSuffix(file) != "trr")
        throw std::runtime_error("ERROR: Called readTRRFile() with wrong extension name.");
    XDRFILE* xd = xdrfile_open(file.c_str(), "r"); // Open a file for xdr I/O ("r" for reading)
    if (xd == NULL)
        throw std::runtime_error("ERROR: Error when opening file: " + file);
    xdrfile_close(xd); // Close the file for xdr I/O
}

void StructureVec3::writeTRRFile(const std::string& file, int step, double time, const std::vector<int>& indices) const {
    if (GetFileSuffix(file) != "trr")
        throw std::runtime_error("ERROR: Called writeTRRFile() with wrong extension name.");
    static int frame = 0; // current frame index in trajectory file
    static std::string file_prev = file; // the filename of previous call
    if (file_prev != file) { // if it is not a same file,
        frame = 0;           // we reset the frame to 0,
        file_prev = file;    // and set current file name as previous one.
    }
    // the number of atoms to save, never change in one trajectory file
    const int nindices = indices.size();
    if (nindices > this->natoms)
        throw std::runtime_error("ERROR: The atom indices in the group to save "
            "trajectory file is out of range.");
    static const int natoms = nindices == 0 ? this->natoms : nindices;
    const std::string mode = frame == 0 ? "w" : "a"; // "w" for writing, "a" for append.
    XDRFILE* xd = xdrfile_open(file.c_str(), mode.c_str()); // Open a file for xdr I/O
    if (xd == NULL)
        throw std::runtime_error("ERROR: Error when opening file to write trr file.");
    // Note: for triclinic box, the boxVectors here in OpenMM reduced form, but
    // no conversion is needed since it can be recognized by VMD and Gromacs.
    static float box[3][3]; // box matrix for trr file
    bool hasPBC = true;
    double volumn = boxVectors[0].dot(boxVectors[1].cross(boxVectors[2]));
    if (volumn == 0.0)
        hasPBC = false; // all box vectors are zero
    else    // for rectangular/cubix only
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                box[i][j] = this->boxVectors[i][j];
    // trr file supports to store positions, velocities, and forces.
    static std::vector<rvec> x(natoms), v(natoms), f(natoms); // rvec is float rvec[3]
    bool hasV = velocities.size() == this->natoms;
    bool hasF = forces.size() == this->natoms;
    if (natoms == this->natoms) { // whole system
        for (int i = 0; i < natoms; ++i)
            for (int j = 0; j < 3; ++j)
                x[i][j] = this->positions[i][j];
        if (hasV)
            for (int i = 0; i < natoms; ++i)
                for (int j = 0; j < 3; ++j)
                    v[i][j] = this->velocities[i][j];
        if (hasF)
            for (int i = 0; i < natoms; ++i)
                for (int j = 0; j < 3; ++j)
                    f[i][j] = this->forces[i][j];
    }
    else {  // the selected atoms
        for (int i = 0; i < natoms; ++i)
            for (int j = 0; j < 3; ++j)
                x[i][j] = this->positions[indices[i]][j];
        if (hasV)
            for (int i = 0; i < natoms; ++i)
                for (int j = 0; j < 3; ++j)
                    v[i][j] = this->velocities[indices[i]][j];
        if (hasF)
            for (int i = 0; i < natoms; ++i)
                for (int j = 0; j < 3; ++j)
                    f[i][j] = this->forces[indices[i]][j];
    }
    // Write or append a frame to xtc file. Return 0 if success.
    // Here, 0.0f is lambda, coupling parameter for free energy methods. (useless here)
    // Here, x.data() return a pointer to the first element, i.e. &x[0].
    int isOK = write_trr(xd, natoms, step, time, 0.0f, (hasPBC ? box : NULL),
        x.data(), (hasV ? v.data() : NULL), (hasF ? f.data() : NULL));
    if (isOK != 0)
        throw std::runtime_error("ERROR: Error when writting data to trr file.");
    xdrfile_close(xd); // Close the file for xdr I/O
    frame++;
}