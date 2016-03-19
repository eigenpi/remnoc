#ifndef _REMNOC_H_
#define _REMNOC_H_

#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <vector>
#include <utility>
#include <functional>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cstring>

using namespace std;

#define TASK_MIGRATION_ENERGY_COST_PER_HOP 1

enum REMAPPING_ALGORITHM { REMAPPING_HEURISTIC = 0, REMAPPING_ANNEALING = 1 };

class GUI_GRAPHICS;

////////////////////////////////////////////////////////////////////////////////
//
// RECT
//
////////////////////////////////////////////////////////////////////////////////

class RECT 
{
 private:
    int _xl, _xr, _yb, _yt;
    int _area;
    double _dist; // from this rect to any other rect in the other set;
 public:
    RECT() { _xl = 0; _xr = 0; _yb = 0; _yt = 0; _area = 0; _dist = 0.; }
    RECT(int xl, int xr, int yb, int yt) : 
        _xl(xl), _xr(xr), _yb(yb), _yt(yt), 
        _area((xr-xl)*(yt-yb)), _dist(0.) {}
    RECT(const RECT &r) : 
        _xl(r._xl), _xr(r._xr), _yb(r._yb), _yt(r._yt),
        _area(r._area), _dist(r._dist) {}
    ~RECT() {}

    int xl() const { return _xl; }
    int xr() const { return _xr; }
    int yb() const { return _yb; }
    int yt() const { return _yt; }
    int area() const { return _area; }
    double dist() const { return _dist; }
    void set_xl(int xl) { _xl = xl; }
    void set_xr(int xr) { _xr = xr; }
    void set_yb(int yb) { _yb = yb; }
    void set_yt(int yt) { _yt = yt; }
    void set_area(int area) { _area = area; }
    void set_dist(double dist) { _dist = dist; }

    void operator=(const RECT &r) {
        _xl = r._xl; _xr = r._xr; _yb = r._yb; _yt = r._yt; 
        _area = r._area;
        _dist = r._dist;
    }
    bool operator==(const RECT &r) const { return (_area == r._area); }
    bool operator!=(const RECT &r) const { return (_area != r._area); }
    bool operator<(const RECT &r)  const { return (_dist < r._dist); }
};

////////////////////////////////////////////////////////////////////////////////
//
// REMNOC
//
////////////////////////////////////////////////////////////////////////////////

// REMapping for Network On Chip; this is the place where everything is
// driven from;

class APPLICATION_NODE;
class APPLICATION_GRAPH;
class RECT;

class REMNOC {

    
 public:
    class TILE {
    private:
        REMNOC *_remnoc; // its owner;
        int _id;
        long _x, _y; // coordinates, address;
        bool _is_broken; // failure of tile recorded here;
        APPLICATION_NODE *_node; // pointer to IP/core node mapped to this tile;
    public:
        TILE(REMNOC *remnoc, long x, long y, long id) :
            _remnoc(remnoc), _x(x), _y(y), _id(id) {
            _is_broken = false;
            _node = 0;
        }
        TILE() : _remnoc(0), _x(0), _y(0), _id(-1) {
            _is_broken = false;
            _node = 0;
        }
        ~TILE() {}
        long x() { return _x; }
        long y() { return _y; }
        int id() { return _id; }
        bool is_broken() { return _is_broken; }
        void set_is_broken(bool val) { _is_broken = val; }
        bool has_node_mapped() { return (_node != 0); }
        APPLICATION_NODE *node() { return _node; }
        void set_id(long id) { _id = id; }
        void set_owner(REMNOC *remnoc) { _remnoc = remnoc; }
        void set(long x, long y, int id) { _x = x; _y = y; _id = id; }
        void set_mapping(APPLICATION_NODE *node) { _node = node; }
        void print_tile();
    };

    
 private:
    // _app_grap will have to be replaced with an array of graphs, when
    // will work with more than one applications concurrently;
    APPLICATION_GRAPH *_app_graph;
    vector<TILE> _tiles;
    // number of columns of tiles/routers of the NOC; for example if the NOC
    // is 4x4 tiles, then _nx = _ny = 4;
    long _nx;
    long _ny; // by default is same as _nx;
    long _tiles_count;
    // _sketch_mapped_tiles has a number of elements equal to the number of,
    // cores and stores id's of tiles which are the new locations of good tiles
    // to which the application cores have to be shuffled within; but this
    // storage is not in order of cores id's; this array represents a polygon;
    // in fact it after remapping it may contain more PEs than the number
    // of nodes, depending on what method I use;
    vector<long> _sketch_mapped_tiles_new; // new;
    // sketch arrays and things used for first main step: finding the new shape;
    vector<int> _sketch_x_coord;
    vector<int> _sketch_y_coord;
    // sets of RECTs used as nodes of the graph which will be partitioned
    // incrementally for area/region migration; these sets will not contain
    // rects with broken tiles (red rects);
    vector<RECT> _rects_black; // regions part of the mapping-shape;
    vector<RECT> _rects_white; // empty regions, available;
    double _x_mass_center_black;
    double _y_mass_center_black;    
    int _displaced_area; // area from initial mapping with new failures;
    int _black_area;
    int _white_area;
    // area of application; minimum needed for correct mapping; now it's the same
    // as _nodes_count basically;
    int _app_area;
    // command line option to record what the user wants to run;
    REMAPPING_ALGORITHM _remapping_algo;
    bool _use_gui;
    // this option should be used with "-use_gui"; if set true by user, 
    // then each user will have to hit "Proceed" button to advance the
    // simulation after every printing interval; default is false;
    int _failures_count;
    int _rng_seed; // seed for random number generator; else is set to 1;
    // storage used during results collection;
    vector<double> _sketch_avg_comm_volume_change;
    vector<double> _sketch_avg_migration_distance;  

 public:
    GUI_GRAPHICS *_gui;
    
 public:
    REMNOC( long nx, long ny);
    REMNOC(); // _nx and _ny will be read from the input mapping_file;
    ~REMNOC() {}
    
    APPLICATION_GRAPH *app_graph() { return _app_graph; };
    GUI_GRAPHICS *gui() { return _gui; };
    void set_gui(GUI_GRAPHICS *gui) { _gui = gui; };
    bool use_gui() { return _use_gui; }
    REMAPPING_ALGORITHM remapping_algo() const { return _remapping_algo; }

    vector<TILE> &tiles() { return _tiles; }
    const vector<TILE> &tiles() const { return _tiles; }
    TILE *tile(long x, long y) {
        assert( x >= 0 && x < _nx && y > 0 && y < _ny);
        long id = x + _nx * y;
        return ( &_tiles[ id]);
    }
    const TILE *tile(long x, long y) const {
        assert( x >= 0 && x < _nx && y > 0 && y < _ny);
        long id = x + _nx * y;
        return ( &_tiles[ id]);
    }
    TILE *tile(long id) {
        assert( id >= 0 && id < _tiles_count);
        return ( &_tiles[ id]);
    }
    long tiles_count() { return _tiles_count; }
    long nx() { return _nx; }
    long ny() { return _ny; }
    void set_nx(long nx) { _nx = nx; }
    void set_ny(long ny) { _ny = ny; }
    void set_mapping_of_tile( int tile_id, APPLICATION_NODE *app_node) {
        _tiles[ tile_id].set_mapping( app_node);
    }
    bool create_application_graphs( int argc, char *argv[]);
    void read_in_mappings( int argc, char *argv[]);
    void build_tiles(); // called by read_in_mappings();
    void clean_tiles_mappings(); // sets tiles' node ptr to zero;
    void inject_failures();
    // SA-based;
    bool run_simulated_annealing_remapping();
    bool run_remapping_single_failure_multiple_times_SA_based();
    // proposed algo;
    bool run_remapping();
    bool run_remapping_single_failure_multiple_times();
    bool run_remapping_two_sequential_failures_multiple_times();
    bool run_remapping_two_simultaneous_failures_multiple_times();
    bool run_remapping_three_simultaneous_failures_multiple_times();
    void run_computation_of_new_map_polygon_1();
    void run_computation_of_new_map_polygon_2();
    void run_assignment_of_cores_to_new_locations();
    void check_remapping();
    long compute_total_migration();
    double compute_total_comm_volume();
    double compute_total_comm_volume_for_annealer(int move_type); // "0" move, "1" swap;
    void build_sketch_xy_coord();
    void build_rects_arrays(); // 1
    void build_rects_arrays_with_unit_area(); // 2
    void compute_mass_center_black(); // update mass center of black rects;
    void compute_mass_center_black_weighted();
    void compute_distances_white_to_black();
    void compute_Euclidean_distances_white_to_mass_center_black();
    int my_white_black_distance( int j, int i);
    double my_Euclidean_white_mass_center_black_distance( int j);
    
    void print_network_tiles();
    void print_sketch_arrays();
    bool parse_command_arguments( int argc, char *argv[]);
    void print_initial_stats( int argc, char *argv[]);
    static bool compare_by_area(const RECT &r1, const RECT &r2) {
        return ( r1.area() > r2.area());
    }
    static bool compare_by_dist(const RECT &r1, const RECT &r2) {
        // want smallest distance, largest area;
        if ( r1.dist() > r2.dist()) {
            return true;
        }
        else if ( (r1.dist() == r2.dist())&&(r1.area() < r2.area()) ) {
            return true;
        }       
        return false;
    }
    static int my_rand_int(int max_val) 
    {
        // return a random integer between [0..max_val];
        // r is a random float value in the range [0,1);
        double r = (double)rand() / ((double)(RAND_MAX)+(double)(1));
        r = r * (max_val + 1);
        return int( r);
    }
};

////////////////////////////////////////////////////////////////////////////////
//
// APPLICATION_NODE
//
////////////////////////////////////////////////////////////////////////////////

// node of application characterization graph; represents basically an
// IP/core of the application already mapped on the NOC;

class APPLICATION_NODE {
 private:
    int _id; // id of itself;
    int _init_tile_id; // initial tile id node was mapped to; need for stats;
    int _tile_id; // id of tile to which is mapped; otherwise -1;
    string _name;
    double _io_comm_volume; // sum of comm_volume of all arcs in or out;
    vector<long> _fanin; // id's of nodes with arcs as destination this node;
    vector<long> _fanout;   
        
 public:
    APPLICATION_NODE() { 
        _id = -1; 
        _tile_id = -1; 
        _init_tile_id = -1; _io_comm_volume = 0.;
    }
    APPLICATION_NODE(int id) : _id(id) { 
        _tile_id = -1; 
        _init_tile_id = -1; _io_comm_volume = 0.;
    }
    APPLICATION_NODE(const APPLICATION_NODE &node) :
        _id(node._id), _tile_id(node._tile_id),
        _init_tile_id(node._init_tile_id), _io_comm_volume(node._io_comm_volume) {}
    ~APPLICATION_NODE() {}
    
    int id() const { return _id; }
    int tile_id() const { return _tile_id; }
    int init_tile_id() const { return _init_tile_id; }
    double io_comm_volume() { return _io_comm_volume; }
    void add_to_io_comm_volume(double delta) { _io_comm_volume += delta; }
    string name() const { return _name; }
    vector<long> &fanin() { return _fanin; }
    vector<long> &fanout() { return _fanout; }
    void set( long id, string name) {
        _id = id;
        _name = name;
    }
    void set_mapped_to( long id) { _tile_id = id; }
    void set_init_mapped_to( long id) { _init_tile_id = id; }
    void add_fanout( long id) { // id is index of core this one fanouts to;
        _fanout.push_back( id);
    }
    void add_fanin( long id) {
        _fanin.push_back( id);
    }
};

////////////////////////////////////////////////////////////////////////////////
//
// APPLICATION_ARC
//
////////////////////////////////////////////////////////////////////////////////

class APPLICATION_ARC {
 private:
    int _id;
    long _src_id; // id of core that is source of this arc;
    long _des_id; // id of core that is destination of this arc;
    double _comm_volume;
        
 public:
    APPLICATION_ARC() { _id = -1; }
    APPLICATION_ARC(int id) : _id(id) {}
    APPLICATION_ARC(const APPLICATION_ARC &arc) : _id(arc._id) {}
    ~APPLICATION_ARC() {}
    
    int id() const { return _id; }
    int src_id() const { return _src_id; }
    int des_id() const { return _des_id; }
    double comm_volume() const { return _comm_volume; }
    void set( long src_id, long des_id, double comm_volume) {
        _src_id = src_id;
        _des_id = des_id;
        _comm_volume = comm_volume;
    }
};

////////////////////////////////////////////////////////////////////////////////
//
// APPLICATION_GRAPH
//
////////////////////////////////////////////////////////////////////////////////

class APPLICATION_GRAPH {
 private:
    int _id;
    long _nodes_count;
    long _arcs_count;
    vector<APPLICATION_NODE> _nodes; // IP/cores of the application;
    vector<APPLICATION_ARC> _arcs; // communication arcs of the application;
    double _max_comm_volume; // maximum communication volume among all nodes; 

 public:
    APPLICATION_GRAPH() { _id = -1; _max_comm_volume = 0; }
    APPLICATION_GRAPH( int id) : _id(id) { _max_comm_volume = 0; }
    ~APPLICATION_GRAPH() {}
    
    int id() const { return _id; }
    long nodes_count() { return _nodes_count; }
    long arcs_count() { return _arcs_count; }
    double max_comm_volume() const { return _max_comm_volume; }
    vector<APPLICATION_NODE> &nodes() { return _nodes; }
    APPLICATION_NODE *get_node(long id) {
        //assert(id >= 0 && id < _nodes_count); 
        return &_nodes[ id]; 
    }
    APPLICATION_ARC *get_arc(long id) {
        return &_arcs[ id]; 
    }
    void build_nodes(long cores_count) {
        _nodes_count = cores_count;
        _nodes.resize( cores_count);
        for (long i = 0; i < _nodes_count; i++) { _nodes[i] = 0; }
    }
    void build_arcs(long arcs_count){ 
        _arcs_count = arcs_count;
        _arcs.resize( arcs_count);
        for (long i = 0; i < _arcs_count; i++) { _arcs[i] = 0; }
    }
    void set_core( long core_id, string core_name) {
        assert( core_id >= 0 && core_id < _nodes_count);
        _nodes[ core_id].set( core_id, core_name);
    }
    void set_arc( long arc_id, long src_id, long des_id, double comm_volume) {
        assert( arc_id >= 0 && arc_id < _arcs_count);
        assert( src_id >= 0 && src_id < _nodes_count);
        assert( des_id >= 0 && des_id < _nodes_count);
        _nodes[ src_id].add_fanout( des_id);
        _nodes[ des_id].add_fanin( src_id);
        // set arc too; this could be avoided by storing this info together with 
        // every fanout of every node;
        _arcs[ arc_id].set( src_id, des_id, comm_volume);
    }
    void compute_max_comm_volume() {
        _max_comm_volume = 0.;
        for ( long i = 0; i < _nodes_count; i++) {
            if ( _max_comm_volume < _nodes[i].io_comm_volume()) {
                _max_comm_volume = _nodes[i].io_comm_volume();
            }
        }
    }
    void print_graph();
};

#endif

