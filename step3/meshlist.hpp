#pragma once
//------------------------------------------------------------------------
#include <vector>
#include "variables.hpp"
//------------------------------------------------------------------------
class MeshList{
  private:
    double mesh_size;
    int m;
    int number_of_mesh;
    std::vector<int> count;
    std::vector<int> index;
    std::vector<int> sorted_buffer;
  public:
    MeshList(void);
    void make_mesh(Variables *vars);
    void set_number_of_atoms(int pn){
      sorted_buffer.resize(pn);
    }
};
//------------------------------------------------------------------------
