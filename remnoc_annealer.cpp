#include "remnoc_annealer.h"
#include <stdio.h>
#include <math.h>

using namespace std;

#define ANNEALER_DELTA 0.01
#define ANNEALER_STOP 0.001
#define ANNEALER_CONGESTION_PENALTY 1000000.0

////////////////////////////////////////////////////////////////////////////////
//
// ANNEALER
//
////////////////////////////////////////////////////////////////////////////////

ANNEALER::ANNEALER( REMNOC *remnoc) : _remnoc(remnoc)
{
    _app_graph = _remnoc->app_graph();
    _nx = _remnoc->nx();
    _ny = _remnoc->ny();    
    _nodes_count_1 = _remnoc->app_graph()->nodes_count() - 1;
    _tiles_count_1 = _remnoc->tiles_count() - 1;
}

void ANNEALER::run_initial_random_assignment()
{
    // assign randomly all nodes to good tiles of the network;
    // () clean up tiles;
    long tiles_count = _tiles_count_1 + 1;
    for ( long j = 0; j < tiles_count; j++) {
        _remnoc->set_mapping_of_tile( j, 0);
    }
    // () random assignment of nodes to good tiles;
    long nodes_count = _nodes_count_1 + 1;
    for ( long i = 0; i < nodes_count; i++) {
        _node_z = _app_graph->get_node( i);
        int t_id = _remnoc->my_rand_int( _tiles_count_1);
        while ( _remnoc->tile( t_id)->is_broken() || // need a good tile;
            _remnoc->tile( t_id)->has_node_mapped()) { // need an empty tile;
            t_id = _remnoc->my_rand_int( _tiles_count_1);
        }

        _node_z->set_mapped_to( t_id);
        _remnoc->set_mapping_of_tile( t_id, _node_z);
    }   
}

////////////////////////////////////////////////////////////////////////////////
//
// a move is moving a cell;
//
////////////////////////////////////////////////////////////////////////////////

double ANNEALER::cost_of_move_node( APPLICATION_NODE *app_node)
{
    // innitially the argument was used to calculate only the cost
    // cost (let's say as sum of bounding boxes of all modes touching this
    // app_node) contribution 'of this app_node";
    // now I am simply computing the total cost;

    double this_cost = 0;

    // () compute total cost or cost-contribution of this node only for
    // faster runtime;
    this_cost = _remnoc->compute_total_comm_volume_for_annealer( 0); // move;

    // () bias a lot the cost if move is to congested clb; this is under the
    // assumption that I work by allowing two nodes to occupy the same tile
    // during the annealing process; I do not do it as I would need
    // to modify the TILE class and could get trapped in local minima;
    // if tile congetsted: cost += ANNEALER_CONGESTION_PENALTY;
    return this_cost;
}

bool ANNEALER::move_one_node( double temperature)
{ 
    // pick up a random application node, generate randomly
    // new PE/tile location, see if it's worth moving it;
    double oldcost = 0, newcost = 0;

    // () select random node;
    int node_i = _remnoc->my_rand_int( _nodes_count_1);
    _node_z = _app_graph->get_node( node_i);
    _node_z_tile_id_1 = _node_z->tile_id();

    // () store oldcost;
    oldcost = cost_of_move_node( _node_z);

    // () generate randomly a new good x,y location (i.e., pickup a random PE/tile);
    _node_z_tile_id_2 = _remnoc->my_rand_int( _tiles_count_1);
    while ( _node_z_tile_id_2 == _node_z_tile_id_1 || // tile 2 shall be diff from current;
        _remnoc->tile( _node_z_tile_id_2)->is_broken() || // need a good tile;
        _remnoc->tile( _node_z_tile_id_2)->has_node_mapped()) { // need an empty tile;

        _node_z_tile_id_2 = _remnoc->my_rand_int( _tiles_count_1);
    }

    // () free initial tile of this node; assign this node to the
    // new location;
    _remnoc->set_mapping_of_tile( _node_z_tile_id_1, 0); // free the PE/tile;
    _node_z->set_mapped_to( _node_z_tile_id_2); // node knows of its new location;
    _remnoc->set_mapping_of_tile( _node_z_tile_id_2, _node_z); // claim the new tile;

    // () after node is moved to new location, compute the new cost;
    newcost = cost_of_move_node( _node_z);

    // () at this time the node has been moved; if it's gonna be decided that
    // the move is not gonna be kept, the caller of this function will have to
    // move back the node - stored in _node_z - to its initial location;
    if ( newcost <= oldcost) {
        return _fix;
    }
    double prob = exp( -(newcost - oldcost) / temperature);
    double val = (rand() / (RAND_MAX + 1.));
    return ( prob > val);
}

void ANNEALER::run_annealer_in_single_node_moves_mode( 
    int seed, double coolio_speed, long moves, double current)
{
    // look first for initial temperature and then do annealing;
    // use only moves that move single nodes;

    srand( seed);
    _stop = _fix = false;
    _range_step = 0.995;
    _range = 2.99;
    long accepted = 0;

    // (1) search for the starting temperature;
    double hot_temperature = 1.;
    while ( fabs( double(accepted)/double(moves) - current) > ANNEALER_DELTA) {
        //printf("\nlooking for HOT temperature = %.2f", hot_temperature);
        //printf("\naccepted %d  %.2f - %.2f = %.2f > %.2f", accepted,
        //  double(accepted)/double(moves), current,
        //  fabs( double(accepted)/double(moves) - current), ANNEALER_DELTA);
        accepted = 0; // reset;
        for ( long m_i = 0; m_i < moves; m_i++) {
            if ( move_one_node( hot_temperature)) {
                accepted ++;
                // move back the node - stored in _node_z - to its initial 
                // location, from where it was relocated by anneal_move_one_node();
                _remnoc->set_mapping_of_tile( _node_z_tile_id_2, 0);
                _node_z->set_mapped_to( _node_z_tile_id_1);
                _remnoc->set_mapping_of_tile( _node_z_tile_id_1, _node_z);
            }
        }
        hot_temperature *= (double(accepted)/double(moves) < current) ? 1.2 : 0.7;
    }
    //printf("\n HOT=%.2f range=%.2f moves=%d ",hot_temperature,_range,moves);

    // (2) do the annealing;
    _fix = true;
    accepted = moves;
    while ( ( double(accepted)/double(moves) > ANNEALER_STOP) && 
            ( _range >= 1.) && // _range is one control;
            ! _stop) {
        accepted = 0;
        for ( long m_i = 0; m_i < moves; m_i++) {
            if ( move_one_node( hot_temperature)) {
                accepted ++;
            } else {
                // if the move was not accepted, then
                // move back the node - stored in _node_z - to its initial 
                // location, from where it was relocated by anneal_move_one_node();
                _remnoc->set_mapping_of_tile( _node_z_tile_id_2, 0);
                _node_z->set_mapped_to( _node_z_tile_id_1);
                _remnoc->set_mapping_of_tile( _node_z_tile_id_1, _node_z);
            }
        }
        _range = _range * _range_step;
        hot_temperature *= coolio_speed; // coolio_speed is another control;
        if ( accepted == 0) { 
            _stop = true; // another control;
        }
        //printf("\n accepted: %d", accepted);
    }
}

////////////////////////////////////////////////////////////////////////////////
//
// a move is swapping two cells;
//
////////////////////////////////////////////////////////////////////////////////

double ANNEALER::cost_of_swap_nodes()
{
    // innitially the argument was used to calculate only the cost
    // cost (let's say as sum of bounding boxes of all modes touching this
    // app_node) contribution 'of this app_node";
    // now I am simply computing the total cost;

    double this_cost = 0;

    // () compute total cost or cost-contribution of this node only for
    // faster runtime;
    this_cost = _remnoc->compute_total_comm_volume_for_annealer( 1); // swap;

    // () bias a lot the cost if move is to congested clb; this is under the
    // assumption that I work by allowing two nodes to occupy the same tile
    // during the annealing process; I do not do it as I would need
    // to modify the TILE class and could get trapped in local minima;
    // if tile congetsted: cost += ANNEALER_CONGESTION_PENALTY;
    return this_cost;
}

bool ANNEALER::move_two_nodes( double temperature)
{ 
    // pick up two random application nodes;
    // see if it's worth swapping them;
    double oldcost = 0, newcost = 0;

    // () select 1st and 2nd random nodes, as _node_u and _node_v;
    int u_i = _remnoc->my_rand_int( _nodes_count_1);
    _node_u = _app_graph->get_node( u_i);
    _node_u_tile_id = _node_u->tile_id();
    int v_i = _remnoc->my_rand_int( _nodes_count_1);
    while ( v_i == u_i) {
        v_i = _remnoc->my_rand_int( _nodes_count_1);
    }
    _node_v = _app_graph->get_node( v_i);
    _node_v_tile_id = _node_v->tile_id();

    // () store oldcost;
    oldcost = cost_of_swap_nodes();

    // () swap the two nodes;
    _remnoc->set_mapping_of_tile( _node_u_tile_id, _node_v);
    _node_v->set_mapped_to( _node_u_tile_id);
    _remnoc->set_mapping_of_tile( _node_v_tile_id, _node_u);
    _node_u->set_mapped_to( _node_v_tile_id);

    // () after nodes are exchanged, compute the new cost;
    newcost = cost_of_swap_nodes();

    // () if it's gonna be decided that the move is not gonna be kept, the caller of 
    // this function will have to reverse the switch;
    if ( newcost <= oldcost) {
        return _fix;
    }
    double prob = exp( -(newcost - oldcost) / temperature);
    double val = (rand() / (RAND_MAX + 1.));
    return ( prob > val);
}

void ANNEALER::run_annealer_in_switch_two_cells_moves_mode( 
    int seed, double coolio_speed, long moves, double current)
{
    // look first for initial temperature and then do annealing;
    // use only moves that move single nodes;
    // this is supposed to be used as a post processing 2nd annealing 
    // process, which should swap nodes within the new mapping region 
    // only, after the 1st annealing process cooled down; in this way
    // total communication cost of the application will be improved even 
    // more but at the expense of extra migration of tasks from initial
    // mapping;

    srand( seed);
    _stop = _fix = false;
    _range_step = 0.995;
    _range = 2.99;
    long accepted = 0;

    // (1) search for the starting temperature;
    double hot_temperature = 1.;
    while ( fabs( double(accepted)/double(moves) - current) > ANNEALER_DELTA) {
        accepted = 0; // reset;
        for ( long m_i = 0; m_i < moves; m_i++) {
            if ( move_two_nodes( hot_temperature)) {
                accepted ++;
                // move back the nodes to their initial locations; reverse 
                // the swap;
                _remnoc->set_mapping_of_tile( _node_u_tile_id, _node_u);
                _node_u->set_mapped_to( _node_u_tile_id);
                _remnoc->set_mapping_of_tile( _node_v_tile_id, _node_v);
                _node_v->set_mapped_to( _node_v_tile_id);
            }
        }
        hot_temperature *= (double(accepted)/double(moves) < current) ? 1.2 : 0.7;
    }

    // (2) do the annealing;
    _fix = true;
    accepted = moves;
    while ( ( double(accepted)/double(moves) > ANNEALER_STOP) && 
            ( _range >= 1.) && // _range is one control;
            ! _stop) {
        accepted = 0;
        for ( long m_i = 0; m_i < moves; m_i++) {
            if ( move_two_nodes( hot_temperature)) {
                accepted ++;
            } else {
                // if the move was not accepted, then
                // move back the nodes to their initial locations; reverse 
                // the swap;
                _remnoc->set_mapping_of_tile( _node_u_tile_id, _node_u);
                _node_u->set_mapped_to( _node_u_tile_id);
                _remnoc->set_mapping_of_tile( _node_v_tile_id, _node_v);
                _node_v->set_mapped_to( _node_v_tile_id);
            }
        }
        _range = _range * _range_step;
        hot_temperature *= coolio_speed; // coolio_speed is another control;
        if ( accepted == 0) { 
            _stop = true; // another control;
        }
    }
}
