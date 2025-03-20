/* -------------------------------------------------------------------------- *
 *                                  QCDyn                                     *
 * -------------------------------------------------------------------------- *
 * This is part of the Quantum Classical Dynamics (QCDyn) program.            *
 *                                                                            *
 * Author: Zhubin Hu @Sun Group @NYU-SH                                       *
 * Last updated: Nov. 23, 2021                                                *
 * -------------------------------------------------------------------------- */

#include "Topology.h"

void Topology::loadTopologyFile(const std::string& file) {
    std::ifstream stream(file.c_str());
    if (!stream)
        throw std::runtime_error("ERROR: Can't open topology file: " + file);
    const std::string suffix = GetFileSuffix(file);
    if (suffix == "top") 
        loadGromacsTopology(stream);
    else if (suffix == "prmtop")
        throw std::runtime_error("ERROR: The Amber topology file is not yet supported.");
    else
        throw std::runtime_error("ERROR: Unkown topology file format: " + suffix);
    stream.close();
}

void Topology::loadGromacsTopology(std::ifstream& topfile) {
    // The following sections must be provided in topology file. And the sections
    // that belong to MoleculeTypes, such as [ atoms ], [ bonds ], [ angles ],
    // [ pairs ], ..., must be after its [ moleculetype ].
    bool defaultsTag = false, atomstypeTag = false, moleculetypeTag = false;
    bool atomsTag = false, moleculesTag = false, systemTag = false;
    // * Start to process the topology file line by line.
    int lineNumber = 0;
    for (std::string line; getline(topfile, line); ) {
        lineNumber++;
        // Lines starting "#include" to include other force filed files are not supported.define specific parameters, e.g.,
        // "#define", or include other force filed files, e.g., "#include".
        if (line.substr(0, 8) == "#include" || line.substr(0, 7) == "#define")
            throw std::runtime_error("ERROR: Sorry, only a standalone topology file is supported now. The #include statements cannot be processed. ");
        // Lines starting "#define" to define specific parameters are not supported.
        // So these lines will be igorned.
        if (line.substr(0, 7) == "#define") continue;
        // Lines starting with semicolon (;) are comments in topology file.
        // A line starting with (;) and (#), or empty, will be ignored.
        // And a line only having space, \t, will be ignored, too.
        // But, you can add comment (start with ;) at the end of a line.
        if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) continue;
        START: // A label for goto, means a new cycle for a new molecule type.
        // If begin with '[', will check it is a valid section or not.
        if (line[0] == '[' && line[line.size() -1] != ']')
            throw std::runtime_error("ERROR: Illegal section line: " + std::to_string(lineNumber) +
                ". Should be surrounded with [] and without any space or character before [ or after ], e.g., [ defaults ]");     
        // * Start to process each section in topology File.
        // Load parameters into ForceField from [ defaults ],
        // which should be the first section as Gromacs request.
        if (line == "[ defaults ]") {
            defaultsTag = true;
            for ( ; getline(topfile, line); ) {
                lineNumber++;
                if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) continue;
                else if (line[0] == '[') break;
                else {
                    std::vector<std::string> fields;
                    SplitLine(fields, line);
                    if (fields.size() < 5)
                        throw std::runtime_error("ERROR: Too few parameters for [ defaults ] in line " + std::to_string(lineNumber));
                    forceFiled.name = "AMBER";
                    forceFiled.vdWType = std::stoi(fields[0]);
                    forceFiled.combinationRule = std::stoi(fields[1]);
                    forceFiled.genPairs = (fields[2] == "yes" ? true : false);
                    forceFiled.LennardJones14Scale = std::stod(fields[3]);
                    forceFiled.Coulomb14Scale = std::stod(fields[4]);
                    if (forceFiled.vdWType != 1)
                        throw std::runtime_error("ERROR: Sorry, only the LJ potential (nbfunc = 1) is supported now.");
                    if (forceFiled.combinationRule != 2)
                        throw std::runtime_error("ERROR: Sorry, only the Lorentz-Berthelot (comb-rule = 2) is supported now.");
                }
            }
        }
        // Load mass and LJ parameters into AtomTypes from [ atomtypes ],
        // which should be before any [ moleculetype ] as Gromacs request.
        if (line == "[ atomtypes ]") {
            atomstypeTag = true;
            for ( ; getline(topfile, line); ) {
                lineNumber++;
                if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) continue;
                else if (line[0] == '[') break;
                else {
                    std::vector<std::string> fields;
                    SplitLine(fields, line);
                    if (fields.size() < 7)
                        throw std::runtime_error("ERROR: Too few parameters for [ atomtypes ] in line " + std::to_string(lineNumber));
                    // Parameters are atomType, atomicNumber, mass, sigma (in nm), epsilon (in kj/mol).
                    atomTypes.push_back({ fields[0], std::stoi(fields[1]), std::stod(fields[2]), std::stod(fields[5]), std::stod(fields[6]) });
                }
            }
        }
        // Load molecule name into MoleculeTypes from [ moleculetype ],
        // which should be before any [ atoms ], [ bonds ], [ angles ], .....
        if (line == "[ moleculetype ]") {
            moleculetypeTag = true;
            for ( ; getline(topfile, line); ) {
                lineNumber++;
                if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) continue;
                else if (line[0] == '[') break;
                else {
                    std::vector<std::string> fields;
                    SplitLine(fields, line);
                    if (fields.size() < 2)
                        throw std::runtime_error("ERROR: Too few parameters for [ moleculetype ] in line " + std::to_string(lineNumber));
                    // We need molecule name only.
                    moleculeTypes.push_back({ fields[0] });
                }
            }
        }
        // The following [ atoms ], [ bonds ], [ angles ], [ dihedrals ], and so
        // on, are belong to current [ moleculetype ] that before them.
        // Load charge of atom into Atoms from [ atoms ].
        if (line == "[ atoms ]") {
            atomsTag = true;
            if (!moleculetypeTag)
                throw std::runtime_error("ERROR: [ atoms ] has to be after [ moleculetype ].");
            for ( ; getline(topfile, line); ) {
                lineNumber++;
                if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) continue;
                else if (line[0] == '[') break;
                else {
                    std::vector<std::string> fields;
                    SplitLine(fields, line);
                    if (fields.size() < 7)
                        throw std::runtime_error("ERROR: Too few parameters for [ moleculetype ] in line " + std::to_string(lineNumber));
                    // Parameters are atomTypeIndex, atomType, and charge.
                    // Here, atomTypeIndex = 0, and will be processed later.
                    moleculeTypes.back().atoms.push_back({ 0, fields[1], std::stod(fields[6]) });
                }
            }
        }
        // Load Harmonic bond parameters into Bonds from [ bonds ].
        if (line == "[ bonds ]") {
            if (!moleculetypeTag)
                throw std::runtime_error("ERROR: [ bonds ] has to be after [ moleculetype ].");
            for ( ; getline(topfile, line); ) {
                lineNumber++;
                if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) continue;
                else if (line[0] == '[') break;
                else {
                    std::vector<std::string> fields;
                    SplitLine(fields, line);
                    if (fields.size() < 5)
                        throw std::runtime_error("ERROR: Too few parameters for [ bonds ] in line " + std::to_string(lineNumber));
                    // Parameters are atomIndex[2], HydrogenTag, distance (in nm), and forceConstant (in kj/mol/nm^2).
                    // Here, HydrogenTag = 0, and will be labled later.
                    moleculeTypes.back().bonds.push_back({ std::stoi(fields[0]), std::stoi(fields[1]), 0, std::stod(fields[3]), std::stod(fields[4]) });
                }
            }
        }
        // Load Harmonic angle parameters into Angles from [ Angles ].
        if (line == "[ angles ]") {
            if (!moleculetypeTag)
                throw std::runtime_error("ERROR: [ angles ] has to be after [ moleculetype ].");
            for ( ; getline(topfile, line); ) {
                lineNumber++;
                if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) continue;
                else if (line[0] == '[') break;
                else {
                    std::vector<std::string> fields;
                    SplitLine(fields, line);
                    if (fields.size() < 6)
                        throw std::runtime_error("ERROR: Too few parameters for [ angles ] in line " + std::to_string(lineNumber));
                    if (std::stoi(fields[3]) != 1)
                        throw std::runtime_error("ERROR: Sorry, only Harmonic angle potential (funct = 1) is supported now.");
                    // Parameters are atomIndex[3], angle (in degree), and forceConstant (in kj/mol/rad^2).
                    // Here, isHAngle = false, and will be labled later.
                    moleculeTypes.back().angles.push_back({ std::stoi(fields[0]), std::stoi(fields[1]),
                        std::stoi(fields[2]), std::stod(fields[4]), std::stod(fields[5]), false });
                }
            }
        }
        // Load dihedral parameters into Dihedrals from [ dihedrals ].
        if (line == "[ dihedrals ]") {
            if (!moleculetypeTag)
                throw std::runtime_error("ERROR: [ dihedrals ] has to be after [ moleculetype ].");
            for ( ; getline(topfile, line); ) {
                lineNumber++;
                if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) continue;
                else if (line[0] == '[') break;
                else {
                    std::vector<std::string> fields;
                    SplitLine(fields, line);
                    if (fields.size() < 8)
                        throw std::runtime_error("ERROR: Too few parameters for [ dihedrals ] in line " + std::to_string(lineNumber));
                    // ! Note: in Gromacs top, the funct = 1 and 9 is different, when
                    // funct = 1, if same dihedral index is defined, only the last one
                    // is valid, while for funct = 9, all definiations of same dihedral
                    // will take effect. However, here, I didn't distinguish them
                    // both funct 1 and 9 dihedrals are loaded (no discard).
                    // since the top file converted from amber prmtop by ParmEd
                    // the funct is always 1 (should be 9). And the OpenMM Python
                    // API also does like this.
                    if (std::stoi(fields[4]) != 1 && std::stoi(fields[4]) != 4 && std::stoi(fields[4]) != 9)
                        throw std::runtime_error("ERROR: Sorry, only Harmonic perodic type dihedral "
                            "potential (funct = 1 (proper) or 4 (improper) or 9 (multiple proper)) are supportted now.");
                    // Parameters are atomIndex[4], periodicity, dihedral (in degree), and forceConstant (in kj/mol/rad^2).
                    moleculeTypes.back().dihedrals.push_back({ std::stoi(fields[0]), std::stoi(fields[1]), std::stoi(fields[2]),
                        std::stoi(fields[3]), std::stoi(fields[7]), std::stod(fields[5]), std::stod(fields[6]) });
                }
            }
        }
        // Load 1-4 interaction into Pairs from [ pairs ].
        // The nonbonded interactions of these atom pairs should be scaled and
        // excluded from normal nonbonded interactions.
        if (line == "[ pairs ]") {
            if (!moleculetypeTag)
                throw std::runtime_error("ERROR: [ pairs ] has to be after [ moleculetype ].");
            for ( ; getline(topfile, line); ) {
                lineNumber++;
                if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) continue;
                else if (line[0] == '[') break;
                else {
                    std::vector<std::string> fields;
                    SplitLine(fields, line);
                    if (fields.size() < 3)
                        throw std::runtime_error("ERROR: Too few parameters for [ dihedrals ] in line " + std::to_string(lineNumber));
                    if (std::stoi(fields[2]) != 1)
                        throw std::runtime_error("ERROR: Sorry, only normal pair interactions (funct = 1) are supported now.");
                    // We only need the atomIndex[2].
                    moleculeTypes.back().pairs.push_back({ std::stoi(fields[0]), std::stoi(fields[1]) });
                }
            }
        }
        // Load parameters for rigid water model into Settles from [ settles ],
        // which defines the first atom of the water molecule. The settle funct
        // is always 1, and the distances between O-H and H-H  must be given.
        // TODO: 4-point water model.
        if (line == "[ settles ]") {
            if (!moleculetypeTag)
                throw std::runtime_error("ERROR: [ bonds ] has to be after [ moleculetype ].");
            for ( ; getline(topfile, line); ) {
                lineNumber++;
                if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) continue;
                else if (line[0] == '[') break;
                else {
                    std::vector<std::string> fields;
                    SplitLine(fields, line);
                    if (fields.size() < 4)
                        throw std::runtime_error("ERROR: Too few parameters for [ settles ] in line " + std::to_string(lineNumber));
                    // Parameters are distances (in nm) between O-H and H-H.
                    moleculeTypes.back().settles.push_back({ std::stod(fields[2]), std::stod(fields[3]) });
                }
            }
        }
        // Load extra exclusions into Exclusions from [ exclusions ].
        // Each line should start with one atom index, followed by one or more
        // atom indices. All non-bonded interactions between the first atom and
        // the other atoms will be excluded.
        if (line == "[ exclusions ]") {
            if (!moleculetypeTag)
                throw std::runtime_error("ERROR: [ exclusions ] has to be after [ moleculetype ].");
            for ( ; getline(topfile, line); ) {
                lineNumber++;
                if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) continue;
                else if (line[0] == '[') break;
                else {
                    std::vector<std::string> fields;
                    SplitLine(fields, line);
                    if (fields.size() < 2)
                        throw std::runtime_error("ERROR: Too few parameters for [ exclusions ] in line " + std::to_string(lineNumber));
                    // Here, the number of atoms is unknown. And comment at the
                    // end of line is allowed. It will stop to read the atom
                    // index when encounter a semicolon (;).
                    std::vector<int> atomIndex;
                    for (auto field : fields) {
                        if (field[0] == ';') break;
                        atomIndex.push_back(std::stoi(field));
                    }
                    moleculeTypes.back().exclusions.push_back({ atomIndex });
                }
            }
        }
        // Load extra constraints for bonds into Constraints from [ constraints ].
        // The format is two atom indices followed by the function type,
        // which can be 1 or 2, and the constraint distance. The only difference
        // between the two types is that type 1 (with chemical bond) is used for
        // generating exclusions and type 2 (without chemical bond) is not.
        // Both types of constraints can be perturbed in free energy calculations
        // by adding a second constraint distance, but is not supported here.
        if (line == "[ constraints ]") {
            if (!moleculetypeTag)
                throw std::runtime_error("ERROR: [ constraints ] has to be after [ moleculetype ].");
            for ( ; getline(topfile, line); ) {
                lineNumber++;
                if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) continue;
                else if (line[0] == '[') break;
                else {
                    std::vector<std::string> fields;
                    SplitLine(fields, line);
                    if (fields.size() < 4)
                        throw std::runtime_error("ERROR: Too few parameters for [ constraints ] in line " + std::to_string(lineNumber));
                    // Parameters are atomIndex[2], useGenExclusion, and distance (in nm).
                    bool useGenExclusion = std::stoi(fields[2]) == 1 ? true : false;
                    moleculeTypes.back().constraints.push_back({ std::stoi(fields[0]), std::stoi(fields[1]), useGenExclusion, std::stod(fields[3]) });
                }
            }
        }
        // Load Harmonic polars parameters into Polars from [ polars ].
        if (line == "[ polars ]") {
            if (!moleculetypeTag)
                throw std::runtime_error("ERROR: [ polars ] has to be after [ moleculetype ].");
            for ( ; getline(topfile, line); ) {
                lineNumber++;
                
                if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) {
                    continue;
                }
                else if (line[0] == '[') { 
                     break;
                }
                else {
                    std::vector<std::string> fields;
                    SplitLine(fields, line);
                    if (fields.size() < 7)
                        throw std::runtime_error("ERROR: Too few parameters for [ moleculetype ] in line " + std::to_string(lineNumber));
                    // Parameters are atomTypeIndex, atomType, and charge.
                    // Here, atomTypeIndex = 0, and will be processed later.
                    moleculeTypes.back().polars.push_back({ 0, fields[1], std::stod(fields[6]) });
                }
            }
        }
        // Load the new vdw model parameters into NewModel from [ newvdw ].
        if (line == "[ newvdw ]") {
            if (!moleculetypeTag)
                throw std::runtime_error("ERROR: [ newvdw ] has to be after [ moleculetype ].");
            for ( ; getline(topfile, line); ) {
                lineNumber++;
                if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) continue;
                else if (line[0] == '[') break;
                else {
                    std::vector<std::string> fields;
                    SplitLine(fields, line);
                    if (fields.size() < 10)
                        throw std::runtime_error("ERROR: Too few parameters for [ moleculetype ] in line " + std::to_string(lineNumber));
                    // Parameters are atomTypeIndex, atomType, and charge.
                    // Here, atomTypeIndex = 0, and will be processed later.
                    moleculeTypes.back().newvdw.push_back({ std::stod(fields[0]), std::stod(fields[1]), std::stod(fields[2]), std::stod(fields[3]),
                                                std::stod(fields[4]), std::stod(fields[5]), std::stod(fields[6]), std::stod(fields[7]),
                                                std::stod(fields[8]), std::stod(fields[9]) });
                }
            }
        }
        // Load system name from [ system ], which should be after all [ moleculetypes ].
        if (line == "[ system ]") {
            systemTag = true;
            for ( ; getline(topfile, line); ) {
                lineNumber++;
                if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) continue;
                else if (line[0] == '[') break;
                else systemName = line;
            }
        }
        // Load the name and number of molecules from [ molecules ],
        // which should be the last section in topology file.
        // Note that the order of molecules listed here must be same with structure
        // file. And the name of molecules must be same with that in [ moleculetype ].
        if (line == "[ molecules ]") {
            moleculesTag = true;
            for ( ; getline(topfile, line); ) {
                lineNumber++;
                if (line[0] == ';' || line[0] == '#' || IsBlankString(line)) continue;
                else if (line[0] == '[') break;
                else {
                    std::vector<std::string> fields;
                    SplitLine(fields, line);
                    if (fields.size() < 2)
                        throw std::runtime_error("ERROR: Too few parameters for [ molecules ] in line " + std::to_string(lineNumber));
                    // Here, moleculeTypeIndex = 0, and will be labled later.
                    molecules.push_back({ fields[0], std::stoi(fields[1]), 0 });
                }
            }
        }
        // Start to process a new molecule type.
        if (line[0] == '[') goto START;
    }
    topfile.close();
    // * Now the processing of topology file line by line is finished.
    // Then we will do some checking for this data set and post processing so as
    // to apply them to OpenMM System more convenient.
    // (1) Check sections in topology file.
    if (!(atomstypeTag && moleculetypeTag && atomsTag && moleculesTag && defaultsTag && systemTag))
        throw std::runtime_error("ERROR: Two few section in your topology file. The [ defaults ], [ atomtypes ], "
            "[ moleculetype ], [ atoms ], [ system ] and [ molecules ] must be provided.");
    // (2) Check if the names in [ molecules ] can be found in [ moleculetype ] and
    //     add moleculeTypeIndex recording the index of MoleculeTypes to Molecules.
    // Note that the names of molecules must be same. However, the order can be
    // different. And one molecule type can be defined several times in [ molecules ]
    // to follow the order of molecules in structure file.
    for (int i = 0; i < molecules.size(); i++)  {
        bool canFind = false;
        for (int j = 0; j < moleculeTypes.size(); j++)
            if (molecules[i].moleculeName == moleculeTypes[j].moleculeName) {
                canFind = true;
                molecules[i].moleculeTypeIndex = j;
                break;
            }
        if (!canFind)
            throw std::runtime_error("ERROR: Can't find the molecule: " + molecules[i].moleculeName +
                " in [ moleculetype ], but exists in the [ molecules ].");
    }
    // (3) Give a label atomTypeIndex in [ atomtypes ] to [ atoms ].
    for (int i = 0; i < moleculeTypes.size(); i++) {
        int count = 0;
        for (int j = 0; j < moleculeTypes[i].atoms.size(); j++)
            for (int k = 0; k < atomTypes.size(); k++)
                if (moleculeTypes[i].atoms[j].atomType == atomTypes[k].atomType) {
                    count++;
                    moleculeTypes[i].atoms[j].atomTypeIndex = k;
                    break;
                }
        if (count != moleculeTypes[i].atoms.size())
            throw std::runtime_error("ERROR: Unknown atomType in [ atoms ] that doesn't exist in [ atomtypes ].");
    }
    // (4) Add a label of HBond or HAngle in bonds/Angles, used for constraints.
    for (int i = 0; i < moleculeTypes.size(); i++) {
        for (int j = 0; j < moleculeTypes[i].bonds.size(); j++) {
            // Note that the atomIndex in [ bonds ] or [ angles ] is from 1.
            // Get the atom index in bonds to find its atomtypeIndex in atoms,
            // then using the atomic number in atomtypes to check H element.
            MoleculeTypes& mtype = moleculeTypes[i];
            MoleculeTypes::Bonds& bond  = mtype.bonds[j];
            MoleculeTypes::Atoms& atomI = mtype.atoms[bond.atomIndex[0] - 1];
            MoleculeTypes::Atoms& atomJ = mtype.atoms[bond.atomIndex[1] - 1];
            AtomTypes& atypeI = atomTypes[atomI.atomTypeIndex];
            AtomTypes& atypeJ = atomTypes[atomJ.atomTypeIndex];
            if (atypeI.atomicNumber == 1) bond.HydrogenTag = 1;
            if (atypeJ.atomicNumber == 1) bond.HydrogenTag = 2;
        }
        for (int j = 0; j < moleculeTypes[i].angles.size(); j++) {
            MoleculeTypes& mtype = moleculeTypes[i];
            MoleculeTypes::Angles& angle = mtype.angles[j];
            MoleculeTypes::Atoms& atomI = mtype.atoms[angle.atomIndex[0] - 1];
            MoleculeTypes::Atoms& atomJ = mtype.atoms[angle.atomIndex[1] - 1];
            MoleculeTypes::Atoms& atomK = mtype.atoms[angle.atomIndex[2] - 1];
            AtomTypes& atypeI = atomTypes[atomI.atomTypeIndex];
            AtomTypes& atypeJ = atomTypes[atomJ.atomTypeIndex];
            AtomTypes& atypeK = atomTypes[atomK.atomTypeIndex];
            // H-X-H or H-O-X (where X is an arbitrary atom) angles can be constrained.
            bool isHXH = atypeI.atomicNumber == 1 && atypeK.atomicNumber == 1;
            bool isHOX = atypeJ.atomicNumber == 8 && (atypeI.atomicNumber == 1 || atypeK.atomicNumber == 1);
            if (isHXH || isHOX) angle.isHAngle = true;
        }
    }
}

int Topology::computeSystemMass() const {
    double systemMass = 0.0;
    for (int i = 0; i < molecules.size(); i++) {
        double moleculeMass = 0.0;
        const Topology::MoleculeTypes& moleculeType = moleculeTypes[molecules[i].moleculeTypeIndex];
        for (int j = 0; j < moleculeType.atoms.size(); j++)
            moleculeMass += atomTypes[moleculeType.atoms[j].atomTypeIndex].mass;
        systemMass += moleculeMass * (double)molecules[i].moleculeNumber;
    }
    return systemMass;
}

void Topology::printTopologyInfo() const {
    // (1) Report the information of force filed that will be used.
    std::cout << "The " <<  forceFiled.name << " force filed will be used for simulation.\n";
    std::cout << "The scale factors being used for 1-4 Coulomb and Lennard-Jones interactions are " <<
        forceFiled.Coulomb14Scale << " and " << forceFiled.LennardJones14Scale <<
        ", respectively.\n";
    // (2) Report the information of system that loaded form topology file.
    std::cout << "The system name is " << systemName << ", and the molecules in it are as following:\n";
    std::cout << "Name\tNumber\ttAtoms\n";
    for (int i = 0; molecules.size(); i++)
        std::cout << molecules[i].moleculeName << "\t" << molecules[i].moleculeNumber <<
            "\t" << moleculeTypes[molecules[i].moleculeTypeIndex].atoms.size() << "\n";
    std::cout << std::endl;
}
