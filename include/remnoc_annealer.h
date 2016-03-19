#ifndef _REMNOC_ANNEALER_H_
#define _REMNOC_ANNEALER_H_

#include "remnoc.h"

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
// ANNEALER
//
////////////////////////////////////////////////////////////////////////////////

class ANNEALER
{
 private:
    int _nx, _ny; // NOC architecture size;
    int _nodes_count_1; // _nodes_count - 1;
    int _tiles_count_1; // _tiles_count - 1;
    // sketch storage;
    APPLICATION_NODE *_node_u;
    APPLICATION_NODE *_node_v;
    APPLICATION_NODE *_node_z;
    int _node_z_tile_id_1; // stores initial tile id of _node_z before move;
    int _node_z_tile_id_2;
    int _node_u_tile_id;
    int _node_v_tile_id;
    int _new_x;
    int _new_y;
    double _oldcost;
    // control mechanisms;
    double _range, _range_step;
    bool _fix;
    bool _stop;
    
 public:
    REMNOC *_remnoc;
    APPLICATION_GRAPH *_app_graph; // copy of pointer to application graph;

 public:
    ANNEALER( REMNOC *remnoc);
    ~ANNEALER() {}

    void run_initial_random_assignment();
    double cost_of_move_node( APPLICATION_NODE *app_node);
    double cost_of_swap_nodes();
    bool move_one_node( double temperature);
    void run_annealer_in_single_node_moves_mode(int seed, 
        double coolio_speed, long moves, double current);
    bool move_two_nodes( double temperature);
    void run_annealer_in_switch_two_cells_moves_mode(
        int seed, double coolio_speed, long moves, double current);
};

#endif
