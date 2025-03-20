#pragma once
#include <cassert>
#include <iosfwd>
#include <set>
#include <string>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <cmath>  
#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <ctime>
#include <sys/time.h>
#include "Vec3.h"

typedef unsigned int AtomIndex;
typedef std::pair<AtomIndex, AtomIndex> AtomPair;
typedef std::vector<AtomPair>  NeighborList;

// Ridiculous O(n^2) version of neighbor list
// for pedagogical purposes and simplicity
// parameter neighborList is automatically clear()ed before 
// neighbors are added
void  computeNeighborListNaive(NeighborList& neighborList,
                              int nAtoms,
                              const std::vector<Vec3>& atomLocations, 
                              std::vector<std::set<int> >& exclusions,
                              const Vec3* periodicBoxVectors,
                              bool usePeriodic,
                              double maxDistance,
                              double minDistance = 0.0,
                              bool reportSymmetricPairs = false
                            );

// O(n) neighbor list method using voxel hash data structure
// parameter neighborList is automatically clear()ed before 
// neighbors are added
void computeNeighborListVoxelHash(
                              NeighborList& neighborList,
                              int nAtoms,
                              const std::vector<Vec3>& atomLocations,
                              std::vector<std::set<int> >& exclusions,
                              const Vec3* periodicBoxVectors,
                              bool usePeriodic,
                              double maxDistance,
                              double minDistance = 0.0,
                              bool reportSymmetricPairs = false
                            );
void computeNoNeighborList(
                          NeighborList& neighborList,
                          int nAtoms);
  
