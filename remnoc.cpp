#include <string>
#include <math.h>
#include <time.h>
#include <stdlib.h>
#include "limits.h"
#include <algorithm>
#include "remnoc_hungarian.h"
#include "remnoc_gui.h"
#include "remnoc.h"
#include "remnoc_annealer.h"

using namespace std;


////////////////////////////////////////////////////////////////////////////////
//
// REMNOC
//
////////////////////////////////////////////////////////////////////////////////

REMNOC::REMNOC( long nx, long ny) :
    _tiles(),
    _nx(nx),
    _ny(ny)
{
    _app_graph = 0;
    _tiles_count = 0;
    _displaced_area = 0;
    _black_area = 0;
    _white_area = 0;
    _x_mass_center_black = 0.;
    _y_mass_center_black = 0.;  
    _app_area = 0;
    _remapping_algo = REMAPPING_HEURISTIC;
    _use_gui = false;
    _failures_count = 0;
    _rng_seed = 1;
    for ( long x = 0; x < _nx; x++) {
        for ( long y = 0; y < _ny; y++) {       
            _tiles.push_back( TILE( this, x, y, _tiles_count));
            _tiles_count++;
        }
    }
    // reading of the application graph and current mapping will be
    // done from the main function; input files will be provided via
    // main's command line arguments for now;
}

REMNOC::REMNOC() :
    _tiles(),
    _nx(0),
    _ny(0),
    _tiles_count(0),
    _displaced_area(0),
    _black_area(0),
    _white_area(0),
    _app_area(0),
    _remapping_algo(REMAPPING_HEURISTIC),
    _use_gui(false),
    _failures_count(0),
    _rng_seed(1)
{
    // _nx and _ny will be read from the mapping_file from first
    // line, and then build() will be called;
    _app_graph = 0;
    _x_mass_center_black = 0.;
    _y_mass_center_black = 0.;
}

bool REMNOC::parse_command_arguments( int argc, char *argv[])
{
    bool result = true;
    int i = 0;

    // () application_graph and current_mapping files are mandatory;
    if (argc < 3) {
        printf("\nUsage:   remnoc app_graph_file current_mapping_file [Options...]\n");
        printf("Example: remnoc app_graph.app mapping.map -failures_count 3\n");
        printf("Options:\n");
        printf("\t[-seed Int]\n");
        printf("\t[-failures_count Int]\n");
        printf("\t[-use_gui]\n");
        printf("\t[-remapping_algo heuristic | annealing]\n");
        exit(1);
    }
    // these two files will be used inside the next two functions;

    // () here we read in the other Options;
    _remapping_algo = REMAPPING_HEURISTIC; // by default we run the new heuristic;
    _use_gui = false; // by default we do not show any graphics;
    i = 3;
    while ( i < argc) {
        if (strcmp (argv[i],"-remapping_algo") == 0) {
            if (argc <= i+1) {
                printf ("Error:  -remapping_algo option requires a string parameter.\n");
                exit (1);
            } 
            if (strcmp(argv[i+1], "heuristic") == 0) {
                _remapping_algo = REMAPPING_HEURISTIC;
            } 
            else if (strcmp(argv[i+1], "annealing") == 0) {
                _remapping_algo = REMAPPING_ANNEALING;
            } else {
                printf("Error:  -routing_algo must be heuristic or annealing.\n");
                exit (1);
            }
            i += 2;
            continue;
        }
        if (strcmp(argv[i],"-use_gui") == 0) {
            _use_gui = true;
            i += 1; // do not take any parameter value;
            continue;
        }
        if (strcmp(argv[i],"-seed") == 0) {
            _rng_seed = atoi( argv[i+1]);
            // just set it right away :);
            srand ( _rng_seed);
            if ((_rng_seed < 1) || (_rng_seed > INT_MAX)) { 
                printf("Error:  -seed value must be between [1 %d].\n", INT_MAX);
                exit(1); 
            }
            i += 2; 
            continue;
        }
        if (strcmp(argv[i],"-failures_count") == 0) {
            _failures_count = atoi( argv[i+1]);
            if ((_failures_count < 0) || (_failures_count > INT_MAX)) { 
                printf("Error:  -failures_count value must be between [0 %d].\n", INT_MAX);
                exit(1); 
            }
            i += 2; 
            continue;
        }
    }
    return result;      
}

void REMNOC::print_initial_stats( int argc, char *argv[])
{
    printf("app_file:                 %s \n", argv[1]);
    printf("mapping_file:             %s \n", argv[2]);
    printf(" NOC dimension:           %d x %d \n", _nx, _ny);
    printf(" IP/cores count:          %d \n", _app_graph->nodes_count());
    printf(" Arcs count:              %d \n", _app_graph->arcs_count());
    printf(" RNG seed:                %d \n", _rng_seed);
    printf(" Failures to inject:      %d \n", _failures_count);
    printf(" Remapping algorithm:     %s \n", (_remapping_algo == REMAPPING_HEURISTIC) ?
        "heuristic" : "annealing");
}

bool REMNOC::create_application_graphs( int argc, char *argv[])
{
    bool result = true;
    int i = 0;
    // () application_graph and current_mapping files are mandatory;
    if (argc < 3) {
        printf("\nUsage:   remnoc app_graph_file current_mapping_file [Options...]\n");
        printf("Example: remnoc app_graph.app mapping.map -failures_count 3\n");
        printf("Options:\n");
        printf("\t[-seed Int]\n");
        printf("\t[-failures_count Int]\n");
        printf("\t[-use_gui]\n");
        printf("\t[-remapping_algo heuristic | annealing]\n");
        exit(1);
    }
    string app_graph_file_name = argv[1];
    ifstream app_graph_ifstream;
    app_graph_ifstream.open( app_graph_file_name.c_str());
    if ( !app_graph_ifstream) {
        printf("\nError: Cannot open app graph file: %s\n",
            app_graph_file_name.c_str());
        exit(1);
    }

    // the application file contains at least one application graph;
    // the format should be:
    // .application ID1 // first line tells me that an application starts;
    // .cores N // there will be following N lines for all "N" cores;
    // 0 DSP    // "0" is the index of this core among all 0..N-1
    // 1 FPGA   // "DSP', "FPGA" are the names of these cores;
    // ...      // cores can be listed in any order of their index; 
    // .arcs M  // there will be following M lines for each directed arc;
    // 0 1 540  // source is core "0", dest is core "1", comm volume is 540;
    // 1 0 300
    // 2 3 250
    // ...
    // .application ID2 // another application; for now we work with one;
    // Notes: these format can be updated, for example we may want to add
    // for each node/core also its area, etc.

    // (b) I assume that the input files are correct for now; do not do any
    // checking and printouts to educate the user;
    long app_count = 0;
    while ( !app_graph_ifstream.eof()) {
        string app_string;
        app_graph_ifstream >> app_string;
        if ( app_string.compare(".application") == 0) {
            app_count ++;
            long app_id;
            app_graph_ifstream >> app_id; // to get "ID1" from .application ID1;

            // () create the new application object;
            _app_graph = new APPLICATION_GRAPH( app_id);

            // () read in cores;
            long this_app_cores_count;
            app_graph_ifstream >> app_string;
            app_graph_ifstream >> this_app_cores_count; // to get "N" from .cores N;
            // build the nodes of the application, which will then be 
            // populated with their names from following N lines;
            _app_graph->build_nodes( this_app_cores_count);

            for ( long i = 0; i < this_app_cores_count; i ++) {
                long this_core_id;
                app_graph_ifstream >> this_core_id; // get "0" from 0 DSP;
                string this_core_name;
                app_graph_ifstream >> this_core_name;
                // record this core;
                _app_graph->set_core( this_core_id, this_core_name); // id, name;
            }

            // () read in arcs;
            long this_app_arcs_count;
            app_graph_ifstream >> app_string;
            app_graph_ifstream >> this_app_arcs_count; // to get "M" from .arcs M;
            // build the arcs of the application, which will then be 
            // populated with their source, dest pairs from following M lines;
            _app_graph->build_arcs( this_app_arcs_count);
            for ( long j = 0; j < this_app_arcs_count; j ++) {
                long src_id, des_id;
                double comm_volume;
                app_graph_ifstream >> src_id;
                app_graph_ifstream >> des_id;
                app_graph_ifstream >> comm_volume;
                _app_graph->set_arc( j, src_id, des_id, comm_volume); // j is index of arc;
                // add comm volume of this arc to its src and des nodes;
                _app_graph->get_node( src_id)->add_to_io_comm_volume( comm_volume);
                _app_graph->get_node( des_id)->add_to_io_comm_volume( comm_volume);
            }
        } 
    }

    // compute and record max comm volume among all nodes;
    _app_graph->compute_max_comm_volume(); // done once only;
    
    app_graph_ifstream.close();
    // sanity check;
    //_app_graph->print_graph();

    return result;
}

void REMNOC::read_in_mappings( int argc, char *argv[])
{
    // this is normally called after create_application_graphs(), which
    // will do some checking;
    bool result = true;
    int i = 0;
    // (a)
    string curr_mapping_file_name = argv[2];
    ifstream curr_map_ifstream;
    curr_map_ifstream.open( curr_mapping_file_name.c_str());
    if ( !curr_map_ifstream) {
        printf("\nError: Cannot open current mapping file: %s\n",
               curr_mapping_file_name.c_str());
        exit(1);
    }

    // the current mapping file contains the mapping of the application graph,
    // which was readin by create_application_graphs();
    // the format should be:
    // .mapping ID1 // first line tells me whose application id this map is;
    // NX NY        // NOC has NX x NY size (eg 4x4);
    // 0 3          // core id 0 is mapped to tile id 3 (location 2,0)
    // 1 7          // core id 1 is mapped to tile id 7 (location 3,1)
    // ...

    // (b) I assume that the input files are correct for now; do not do any
    // checking and printouts to educate the user;
    string string_t;
    curr_map_ifstream >> string_t;
    long map_id;
    curr_map_ifstream >> map_id; // to get "ID1" from .mapping ID1;
    curr_map_ifstream >> _nx;
    curr_map_ifstream >> _ny;
    _tiles_count = _nx * _ny;

    // () create the tiles for the new network NX x NY object;
    build_tiles();

    // () read in mapping;
    while ( !curr_map_ifstream.eof()) {
        long core_id, tile_id;
        curr_map_ifstream >> core_id;
        curr_map_ifstream >> tile_id;
        assert(core_id >=0 && core_id < _app_graph->nodes_count());
        assert(tile_id >=0 && tile_id < _tiles_count);
        APPLICATION_NODE *app_node = _app_graph->get_node( core_id);
        _tiles[ tile_id].set_mapping( app_node); // tile knows what node;
        app_node->set_mapped_to( tile_id); // node knows to what tile;
        app_node->set_init_mapped_to( tile_id); // record for stats calculation;
    }

    // (c) prepare the sketch arrays for new mapping as well here;
    long nodes_count = _app_graph->nodes_count();
    _sketch_mapped_tiles_new.resize( nodes_count);
    for ( long j = 0; j < nodes_count; j++) {
        _sketch_mapped_tiles_new[ j] = 0; // reset;
    }

    // (d) close mapping file;
    curr_map_ifstream.close();
    // sanity check;
    //print_network_tiles(); // will print tiles together with mapping;
}

void REMNOC::build_sketch_xy_coord()
{
    // call after failure injection (from therein actually) and 
    // before remapping;
    _sketch_x_coord.clear();
    _sketch_y_coord.clear();
    
    vector<bool> x_coord_has_been_recorded;
    vector<bool> y_coord_has_been_recorded;
    x_coord_has_been_recorded.resize( _nx);
    for ( int k = 0; k < _nx; k++) {
        x_coord_has_been_recorded[k] = false;
    }
    y_coord_has_been_recorded.resize( _ny);
    for ( int k = 0; k < _ny; k++) {
        y_coord_has_been_recorded[k] = false;
    }

    // ()
    _sketch_x_coord.push_back( 0);
    _sketch_y_coord.push_back( 0);
    for ( int j = 0; j < _ny; j++) {
        for ( int i = 0; i < _nx; i++) {
            int t_id = i + j * _nx;
            int t_id_east = t_id + 1;
            int t_id_north = t_id + _nx;
            if ( i < _nx - 1) {
                if ( _tiles[t_id].has_node_mapped() != _tiles[t_id_east].has_node_mapped() ||
                    _tiles[t_id].is_broken() != _tiles[t_id_east].is_broken()) {
                    if ( x_coord_has_been_recorded[ i+1] == false) {
                        x_coord_has_been_recorded[ i+1] = true;
                        _sketch_x_coord.push_back( i+1); // record new x coord;
                    }
                }
            }
            if ( j < _ny - 1) {
                if ( _tiles[t_id].has_node_mapped() != _tiles[t_id_north].has_node_mapped() ||
                    _tiles[t_id].is_broken() != _tiles[t_id_north].is_broken()) {
                    if ( y_coord_has_been_recorded[ j+1] == false) {
                        y_coord_has_been_recorded[ j+1] = true;
                        _sketch_y_coord.push_back( j+1); // record new y coord;
                    }
                }
            }
        }
    }
    _sketch_x_coord.push_back( _nx);
    _sketch_y_coord.push_back( _ny);

    // () sort;
    sort( _sketch_x_coord.begin(), _sketch_x_coord.end());
    sort( _sketch_y_coord.begin(), _sketch_y_coord.end());

    // () sanity check;
    //print_sketch_arrays();
}

void REMNOC::build_rects_arrays()
{
    _black_area = 0;
    _white_area = 0;
    
    int x_count = _sketch_x_coord.size();
    int y_count = _sketch_y_coord.size();
    for ( int j = 0; j < y_count - 1; j++) {
        for ( int i = 0; i < x_count - 1; i++) {
            int xl = _sketch_x_coord[ i];
            int xr = _sketch_x_coord[ i+1];
            int yb = _sketch_y_coord[ j];
            int yt = _sketch_y_coord[ j+1];
            int t_id = xl + yb * _nx;
            bool is_red = _tiles[t_id].is_broken();
            bool is_black = ( !is_red && _tiles[t_id].has_node_mapped());
            bool is_white = ( !is_red && !_tiles[t_id].has_node_mapped());
            if ( !is_red) assert( is_black != is_white);
            if ( is_black) {
                _rects_black.push_back( RECT(xl,xr,yb,yt));
                _black_area += (xr - xl) + (yt - yb);
            }
            if ( is_white) {
                _rects_white.push_back( RECT(xl,xr,yb,yt));
                _white_area += (xr - xl) + (yt - yb);
            }
        }
    }
    // () sanity check;
    //print_sketch_arrays();    
}

void REMNOC::build_rects_arrays_with_unit_area()
{
    _black_area = 0;
    _white_area = 0;
    _rects_black.clear();
    _rects_white.clear();   
    for ( int j = 0; j < _ny; j++) {
        for ( int i = 0; i < _nx; i++) {
            int t_id = i + j * _nx;
            bool is_red = _tiles[t_id].is_broken();
            bool is_black = ( !is_red && _tiles[t_id].has_node_mapped());
            bool is_white = ( !is_red && !_tiles[t_id].has_node_mapped());
            if ( !is_red) assert( is_black != is_white);
            if ( is_black) {
                _rects_black.push_back( RECT(i, i+1, j, j+1));
                _black_area ++;
            }
            if ( is_white) {
                _rects_white.push_back( RECT(i, i+1, j, j+1));
                _white_area ++;
            }
        }
    }
}

void REMNOC::compute_mass_center_black()
{
    // all rects have mass = area = 1;
    _x_mass_center_black = 0.;
    _y_mass_center_black = 0.;
    int rects_black_count = _rects_black.size();
    for ( int i = 0; i < rects_black_count; i ++) {
        _x_mass_center_black += 
            ( double(_rects_black[i].xl() + _rects_black[i].xr()) / 2 );
        _y_mass_center_black += 
            ( double(_rects_black[i].yb() + _rects_black[i].yt()) / 2 );
    }
    _x_mass_center_black /= rects_black_count;
    _y_mass_center_black /= rects_black_count;

    //printf("\n Blacks mass center x_c,y_c = %.1f, %.1f",
    //  _x_mass_center_black, _y_mass_center_black);
}

void REMNOC::compute_mass_center_black_weighted()
{
    // all rects have mass = area = 1;
    _x_mass_center_black = 0.;
    _y_mass_center_black = 0.;
    double total_weight = 0.;
    double upper_bound_volume = _app_graph->max_comm_volume() + 1;
    int rects_black_count = _rects_black.size();
    for ( int i = 0; i < rects_black_count; i ++) {
        int xl = _rects_black[i].xl();
        int yb = _rects_black[i].yb();
        long t_id = xl + yb * _ny;
        APPLICATION_NODE *app_node = _tiles[ t_id].node();
        double r_weight = (app_node != 0) ?
            (app_node->io_comm_volume()) : 1.;
        total_weight += r_weight;
        _x_mass_center_black += 
            ( double(_rects_black[i].xl() + _rects_black[i].xr()) / 2 ) * r_weight;
        _y_mass_center_black += 
            ( double(_rects_black[i].yb() + _rects_black[i].yt()) / 2 ) * r_weight;
    }
    _x_mass_center_black /= total_weight;
    _y_mass_center_black /= total_weight;

    //printf("\n Blacks mass center x_c,y_c = %.1f, %.1f",
    //  _x_mass_center_black, _y_mass_center_black);
}

void REMNOC::print_sketch_arrays()
{
    printf("\nx_coord:");
    for ( int i = 0; i < _sketch_x_coord.size(); i ++) {
        printf(" %d", _sketch_x_coord[i]);
    }
    printf("\ny_coord:");
    for ( int i = 0; i < _sketch_y_coord.size(); i ++) {
        printf(" %d", _sketch_y_coord[i]);
    }
    printf("\nrects_black:");
    for ( int i = 0; i < _rects_black.size(); i ++) {
        printf(" (%d,%d,%d,%d):%d ",
            _rects_black[i].xl(),_rects_black[i].xr(),
            _rects_black[i].yb(),_rects_black[i].yt(), 
            _rects_black[i].area());
    }
    printf("\nrects_white:");
    for ( int i = 0; i < _rects_white.size(); i ++) {
        printf(" (%d,%d,%d,%d):%d dist:%.1f ",
            _rects_white[i].xl(),_rects_white[i].xr(),
            _rects_white[i].yb(),_rects_white[i].yt(),
            _rects_white[i].area(), _rects_white[i].dist());
    }
}

void REMNOC::build_tiles()
{
    // by this time _tiles_count should have been read from mapping file;
    _tiles.resize( _tiles_count);
    int id = 0;
    for ( long j = 0; j < _ny; j++) {
        for ( long i = 0; i < _nx; i++) {
            id = i + j * _nx;
            _tiles[ id].set(i, j, id); // x,y,id;
            _tiles[ id].set_owner( this);
        }
    }
}

void REMNOC::clean_tiles_mappings()
{
    for ( long j = 0; j < _tiles_count; j++) {
        _tiles[ j].set_mapping( 0);
    }
}

void REMNOC::print_network_tiles()
{
    // debug;
    printf("\nTILES of NOC %d x %d", _nx, _ny);
    for ( long j = _ny-1; j >= 0; j--) {
        printf("\n");
        for ( long i = 0; i < _nx; i++) {
            long t_id = i + j * _nx;
            if ( _tiles[t_id].is_broken()) {
                printf("  X");
            } else {
                if ( _tiles[t_id].has_node_mapped()) {
                    printf("%3d", _tiles[t_id].node()->id());
                } else {
                    printf("  -");
                }
            }
        }
    }
    //for ( long i = 0; i < _tiles_count; i ++) {
    //printf("\n%d(%d, %d) %s", 
    //     _tiles[i].id(), _tiles[i].x(), _tiles[i].y(),
    //     (_tiles[i].is_broken() ? "FAIL " : "  OK "));
    //if ( _tiles[i].has_node_mapped()) {
    //  printf(" HAS MAPPED CORE %d < %s >", _tiles[i].node()->id(),
    //      _tiles[i].node()->name().c_str());
    //}
    //}
}

void REMNOC::inject_failures()
{
    // mark as broken/failure selected tiles according to 
    // requested distribution;
    // (a)
    long nodes_count = _app_graph->nodes_count();
    int white_space = _tiles_count - nodes_count;
    if ( white_space < _failures_count) {
        printf("\nWarning:  Required failures count %d larger than available white space %d",
           _failures_count, _white_area);
        printf("\nAdjusting failures_count...");
        _failures_count = white_space;
    }
    //printf("\nInjecting  %d  random failures into  %d  total tiles...",
    //     _failures_count, (_nx*_ny)); 
    
    // (b) injection; a number of _failures_count are marked as broken
    // out of the nodes of the application, NOT out of the total tiles
    // count;
    int counter = 0;
    while ( counter < _failures_count) {
        int t_id = _app_graph->get_node( my_rand_int(int(nodes_count)))->tile_id();
        if ( !_tiles[ t_id].is_broken() && 
            _tiles[ t_id].has_node_mapped()) { // keep trying till break req. number;
            _tiles[ t_id].set_is_broken( true);
            counter ++;
        }
    }


    // (c) compute _displaced_area from current/init mapping by failures;
    _displaced_area = 0;
    for ( long i = 0; i < _tiles_count; i ++) {
        if ( _tiles[i].is_broken()) {
            // TODO: this will have to be done for evey application
            // when I will work with more concurrent apps;
            _displaced_area ++;
        }
    }

    // (d) start populating sketch x,y arrays;
    build_sketch_xy_coord();
}

////////////////////////////////////////////////////////////////////////////////
//
// SA-based algo;
//
////////////////////////////////////////////////////////////////////////////////

bool REMNOC::run_simulated_annealing_remapping()
{
    char msg[BUFFER_SIZE];

    ANNEALER my_annealer( this);

    // () random initial assignment;
    my_annealer.run_initial_random_assignment();

    // () entertain user;
    if ( use_gui()) {
        sprintf( msg, "Initial random mapping for annealer");
        gui()->update_screen( PRIORITY_MAJOR, msg, TILES);
    }   
    
    // () run simulated annealing;
    // 1st post processing annealing: move is moving nodes;
    long nodes_count = _app_graph->nodes_count();
    my_annealer.run_annealer_in_single_node_moves_mode(
        _rng_seed, // seed;
        0.99, // coolio_speed;
        long(2000 * nodes_count), // moves;
        0.02); // current;

    // 2nd post processing annealing: move is swapping nodes;
    //my_annealer.run_annealer_in_switch_two_cells_moves_mode(
    //  _rng_seed, // seed;
    //  0.99, // coolio_speed;
    //  long(2000 * nodes_count), // moves;
    //  0.02); // current;
}

bool REMNOC::run_remapping_single_failure_multiple_times_SA_based()
{
    // this is essentially run_simulated_annealing_remapping() run multiple 
    // times to collect results:
    // 1 - a single failure is injected in every PE
    //     on which the application is initially mapped to;
    // 2 - run SA-based remapping, store result;
    // 3 - restore initial mapping, inject new failure; goto 2;

    char msg[BUFFER_SIZE];
    double comm_volume_init = compute_total_comm_volume();
    _sketch_avg_migration_distance.clear();
    _sketch_avg_comm_volume_change.clear();

    // (1)
    long nodes_count = _app_graph->nodes_count();
    for ( long i = 0; i < nodes_count; i++) {
        //printf(" %d", i);
        // () GUI stuff: entertain user before;
        if ( use_gui()) {
            sprintf( msg, "Initial mapping");
            gui()->update_screen( PRIORITY_MAJOR, msg, TILES);
        }       

        // () inject failure into this node's tile/PE;
        _failures_count = 1;
        long t_id_failure = _app_graph->get_node( i)->tile_id();
        _tiles[  t_id_failure].set_is_broken( true);


        ANNEALER my_annealer( this);

        // () random initial assignment;
        my_annealer.run_initial_random_assignment();

        // () run simulated annealing;
        // 1st post processing annealing: move is moving nodes;
        long nodes_count = _app_graph->nodes_count();
        my_annealer.run_annealer_in_single_node_moves_mode(
            time(NULL), // seed;
            0.99, // coolio_speed;
            long(3000 * nodes_count), // moves;
            0.02); // current;

        // () GUI stuff: entertain user after;
        if ( use_gui()) {
            sprintf( msg, "Mapping after failure injection %d", i);
            gui()->update_screen( PRIORITY_MAJOR, msg, TILES);
        }

        // () compute results and store;
        double total_migration_distance = double(compute_total_migration());
        double comm_volume_final = compute_total_comm_volume();
        _sketch_avg_migration_distance.push_back( total_migration_distance);
        double total_comm_volume_change = 100 * double(comm_volume_final - comm_volume_init)/
            double(comm_volume_init);
        _sketch_avg_comm_volume_change.push_back( total_comm_volume_change);

        // () restore initial mapping to prepare for next single failure injection;
        // _rects_black, _rects_white, _black_area, _white_area are reset
        // by run_computation_of_new_map_polygon_2();
        _sketch_mapped_tiles_new.clear();
        _sketch_mapped_tiles_new.resize( nodes_count);
        for ( long sketch_i = 0; sketch_i < nodes_count; sketch_i++) {
            _sketch_mapped_tiles_new[ sketch_i] = 0;
        }
        for ( long j = 0; j < _tiles_count; j++) {
            _tiles[ j].set_mapping( 0);
            _tiles[ j].set_is_broken( false);
        }   
        for ( long j = 0; j < nodes_count; j++) {
            APPLICATION_NODE *app_node = _app_graph->get_node( j);
            long t_id = app_node->init_tile_id();
            app_node->set_mapped_to( t_id); // node knows to what tile;
            _tiles[ t_id].set_mapping( app_node);
        }
    }


    // (2) entertain user; compute overall results with std_dev.;
    long single_failures_count = _sketch_avg_migration_distance.size();
    double avg_migration_distance = 0;
    double avg_migration_distance_STDDEV = 0;   
    double avg_comm_volume_change = 0;
    double avg_comm_volume_change_STDDEV = 0;   
    for ( long i = 0; i < single_failures_count; i++) {
        avg_migration_distance += _sketch_avg_migration_distance[ i];
        avg_comm_volume_change += _sketch_avg_comm_volume_change[ i];
        printf ("\n %d migr_dist: %.2f comm_vol_change: %.2f", i,
            _sketch_avg_migration_distance[ i], _sketch_avg_comm_volume_change[ i]);
    }
    avg_migration_distance /= single_failures_count;
    avg_comm_volume_change /= single_failures_count;
    for ( long i = 0; i < single_failures_count; i++) {
        avg_migration_distance_STDDEV = avg_migration_distance_STDDEV +
            pow( (_sketch_avg_migration_distance[ i] - avg_migration_distance), 2);
        avg_comm_volume_change_STDDEV = avg_comm_volume_change_STDDEV +
            pow( (_sketch_avg_comm_volume_change[ i] - avg_comm_volume_change), 2);
    }
    avg_migration_distance_STDDEV /= single_failures_count;
    avg_comm_volume_change_STDDEV /= single_failures_count;
    avg_migration_distance_STDDEV = sqrt(avg_migration_distance_STDDEV);
    avg_comm_volume_change_STDDEV = sqrt(avg_comm_volume_change_STDDEV);

    printf ("\nSINGLE FAILURE INJECTION RESULTS");
    printf ("\n Avg. total migration distance: %.2f (STDDEV: %.2f)",
            avg_migration_distance, avg_migration_distance_STDDEV);
    printf ("\n Avg. comm volume change: %.2f (STDDEV: %.2f)\n",
            avg_comm_volume_change, avg_comm_volume_change_STDDEV);
}

////////////////////////////////////////////////////////////////////////////////
//
// proposed algo;
//
////////////////////////////////////////////////////////////////////////////////

bool REMNOC::run_remapping()
{
    // application is remapped so that:
    // - broken tiles/cores are not used
    // - remapping is done on a set of tiles overlapping as much as
    //   possible with the original "map-shape";
    // - new locations of every IP/core should be as close as possible
    //   to their initial locations in the original map: minimize
    //   the total task migration amount for energy;

    // (a) get the new mapping shape/polygon region (ideally convex);
    //run_computation_of_new_map_polygon_1();
    run_computation_of_new_map_polygon_2();

    // (b) assign cores to tiles inside the new mapping "shape";
    run_assignment_of_cores_to_new_locations();
}

bool REMNOC::run_remapping_single_failure_multiple_times()
{
    // this is essentially run_remapping() run multiple times to
    // collect results:
    // 1 - a single failure is injected in every PE
    //     on which the application is initially mapped to;
    // 2 - run remapping, store result;
    // 3 - restore initial mapping, inject new failure; goto 2;

    char msg[BUFFER_SIZE];
    double comm_volume_init = compute_total_comm_volume();
    _sketch_avg_migration_distance.clear();
    _sketch_avg_comm_volume_change.clear();

    // (1)
    long nodes_count = _app_graph->nodes_count();
    for ( long i = 0; i < nodes_count; i++) {
        //printf(" %d", i);
        // () GUI stuff: entertain user before;
        if ( use_gui()) {
            sprintf( msg, "Initial mapping");
            gui()->update_screen( PRIORITY_MAJOR, msg, TILES);
        }       

        // () inject failure into this node's tile/PE;
        _failures_count = 1;
        long t_id_failure = _app_graph->get_node( i)->tile_id();
        _tiles[  t_id_failure].set_is_broken( true);

        // () get the new mapping shape/polygon region (ideally convex);
        //run_computation_of_new_map_polygon_1();
        run_computation_of_new_map_polygon_2();

        // () assign cores to tiles inside the new mapping "shape";
        run_assignment_of_cores_to_new_locations();

        // () GUI stuff: entertain user after;
        if ( use_gui()) {
            sprintf( msg, "Mapping after failure injection %d", i);
            gui()->update_screen( PRIORITY_MAJOR, msg, TILES);
        }

        // () compute results and store;
        double total_migration_distance = double(compute_total_migration());
        double comm_volume_final = compute_total_comm_volume();
        _sketch_avg_migration_distance.push_back( total_migration_distance);
        double total_comm_volume_change = 100 * double(comm_volume_final - comm_volume_init)/
            double(comm_volume_init);
        _sketch_avg_comm_volume_change.push_back( total_comm_volume_change);

        // () restore initial mapping to prepare for next single failure injection;
        // _rects_black, _rects_white, _black_area, _white_area are reset
        // by run_computation_of_new_map_polygon_2();
        _sketch_mapped_tiles_new.clear();
        _sketch_mapped_tiles_new.resize( nodes_count);
        for ( long sketch_i = 0; sketch_i < nodes_count; sketch_i++) {
            _sketch_mapped_tiles_new[ sketch_i] = 0;
        }
        for ( long j = 0; j < _tiles_count; j++) {
            _tiles[ j].set_mapping( 0);
            _tiles[ j].set_is_broken( false);
        }   
        for ( long j = 0; j < nodes_count; j++) {
            APPLICATION_NODE *app_node = _app_graph->get_node( j);
            long t_id = app_node->init_tile_id();
            app_node->set_mapped_to( t_id); // node knows to what tile;
            _tiles[ t_id].set_mapping( app_node);
        }
    }


    // (2) entertain user; compute overall results with std_dev.;
    long single_failures_count = _sketch_avg_migration_distance.size();
    double avg_migration_distance = 0;
    double avg_migration_distance_STDDEV = 0;   
    double avg_comm_volume_change = 0;
    double avg_comm_volume_change_STDDEV = 0;   
    for ( long i = 0; i < single_failures_count; i++) {
        avg_migration_distance += _sketch_avg_migration_distance[ i];
        avg_comm_volume_change += _sketch_avg_comm_volume_change[ i];
        printf ("\n %d migr_dist: %.2f comm_vol_change: %.2f", i,
            _sketch_avg_migration_distance[ i], _sketch_avg_comm_volume_change[ i]);
    }
    avg_migration_distance /= single_failures_count;
    avg_comm_volume_change /= single_failures_count;
    for ( long i = 0; i < single_failures_count; i++) {
        avg_migration_distance_STDDEV = avg_migration_distance_STDDEV +
            pow( (_sketch_avg_migration_distance[ i] - avg_migration_distance), 2);
        avg_comm_volume_change_STDDEV = avg_comm_volume_change_STDDEV +
            pow( (_sketch_avg_comm_volume_change[ i] - avg_comm_volume_change), 2);
    }
    avg_migration_distance_STDDEV /= single_failures_count;
    avg_comm_volume_change_STDDEV /= single_failures_count;
    avg_migration_distance_STDDEV = sqrt(avg_migration_distance_STDDEV);
    avg_comm_volume_change_STDDEV = sqrt(avg_comm_volume_change_STDDEV);

    printf ("\nSINGLE FAILURE INJECTION RESULTS");
    printf ("\n Avg. total migration distance: %.2f (STDDEV: %.2f)",
            avg_migration_distance, avg_migration_distance_STDDEV);
    printf ("\n Avg. comm volume change: %.2f (STDDEV: %.2f)\n",
            avg_comm_volume_change, avg_comm_volume_change_STDDEV);
}

bool REMNOC::run_remapping_two_sequential_failures_multiple_times()
{
    // this is essentially run_remapping() run multiple times to
    // collect results:
    // 1 - first a single failure is injected in every PE
    //     on which the application is initially mapped to;
    // 2 - run remapping;
    // 3 - for evey first single failure we then inject a second
    //     failure in any possible remaining core; compute avg.
    //     for all pairs; store result;
    // 3 - restore initial mapping, inject new first failure in the
    //     next core; goto 2;
    // we collect this way results for two failure sequential injections;
    // runtime will be longer...;

    char msg[BUFFER_SIZE];
    double comm_volume_init = compute_total_comm_volume();
    _sketch_avg_migration_distance.clear();
    _sketch_avg_comm_volume_change.clear();

    // (1)
    long nodes_count = _app_graph->nodes_count();
    for ( long i1 = 0; i1 < nodes_count; i1++) {
        // () GUI stuff: entertain user before;
        if ( use_gui()) {
            sprintf( msg, "Initial mapping");
            gui()->update_screen( PRIORITY_MAJOR, msg, TILES);
        }       

        // () inject first failure into this node's tile/PE;
        _failures_count = 1;
        long t_id_failure_1 = _app_graph->get_node( i1)->tile_id();
        _tiles[ t_id_failure_1].set_is_broken( true);

        // () get the new mapping shape/polygon region (ideally convex);
        //run_computation_of_new_map_polygon_1();
        run_computation_of_new_map_polygon_2();
        // () assign cores to tiles inside the new mapping "shape";
        run_assignment_of_cores_to_new_locations();

        // () GUI stuff: entertain user after;
        if ( use_gui()) {
            sprintf( msg, "Mapping after first failure injection %d", i1);
            gui()->update_screen( PRIORITY_MAJOR, msg, TILES);
        }

        // inject second failure in any other remaining core; this is done
        // with restoration of the first remapping;
        vector <long> sketch_mapped_tiles_1;
        // first store the result of remapping after first injection;
        for ( long k = 0; k < nodes_count; k++) {
            sketch_mapped_tiles_1.push_back( _app_graph->get_node( k)->tile_id());
        }
        double total_migration_distance = 0.;
        double comm_volume_final = 0.;
        double total_comm_volume_change= 0.;
        for ( long i2 = i1 + 1; i2 < nodes_count; i2++) {
        
            // () GUI stuff: entertain user after;
            if ( use_gui()) {
                sprintf( msg, "Mapping before second failure injection %d %d",i1,i2);
                gui()->update_screen( PRIORITY_MAJOR, msg, TILES);
            }

            // inject second failure;
            long t_id_failure_2 = _app_graph->get_node( i2)->tile_id();
            _tiles[ t_id_failure_2].set_is_broken( true);

            // () get the new mapping shape/polygon region (ideally convex);
            //run_computation_of_new_map_polygon_1();
            run_computation_of_new_map_polygon_2();
            // () assign cores to tiles inside the new mapping "shape";
            run_assignment_of_cores_to_new_locations();

            // () GUI stuff: entertain user after;
            if ( use_gui()) {
                sprintf( msg, "Mapping after second failure injection %d %d",i1,i2);
                gui()->update_screen( PRIORITY_MAJOR, msg, TILES);
            }

            // record the result of this i1 i2 failures sequence;
            double total_migration_distance = double(compute_total_migration());
            double comm_volume_final = compute_total_comm_volume();
            double total_comm_volume_change = 100 * double(comm_volume_final - comm_volume_init)/
                double(comm_volume_init);
            _sketch_avg_migration_distance.push_back( total_migration_distance);
            _sketch_avg_comm_volume_change.push_back( total_comm_volume_change);

            // restore first remapping for i1;
            _sketch_mapped_tiles_new.clear();
            _sketch_mapped_tiles_new.resize( nodes_count);
            for ( long sketch_i = 0; sketch_i < nodes_count; sketch_i++) {
                _sketch_mapped_tiles_new[ sketch_i] = 0;
            }
            for ( long j = 0; j < _tiles_count; j++) {
                _tiles[ j].set_mapping( 0);
            }
            // note I am backing up only i2 failure;
            _tiles[ t_id_failure_2].set_is_broken( false);  
            for ( long j = 0; j < nodes_count; j++) {
                APPLICATION_NODE *app_node = _app_graph->get_node( j);
                long t_id = sketch_mapped_tiles_1[ j];
                app_node->set_mapped_to( t_id); // node knows to what tile;
                _tiles[ t_id].set_mapping( app_node);
            }
        }


        // () restore initial mapping to prepare for next i1 first failure injection;
        // _rects_black, _rects_white, _black_area, _white_area are reset
        // by run_computation_of_new_map_polygon_2();
        _sketch_mapped_tiles_new.clear();
        _sketch_mapped_tiles_new.resize( nodes_count);
        for ( long sketch_i = 0; sketch_i < nodes_count; sketch_i++) {
            _sketch_mapped_tiles_new[ sketch_i] = 0;
        }
        for ( long j = 0; j < _tiles_count; j++) {
            _tiles[ j].set_mapping( 0);
            _tiles[ j].set_is_broken( false);
        }   
        for ( long j = 0; j < nodes_count; j++) {
            APPLICATION_NODE *app_node = _app_graph->get_node( j);
            long t_id = app_node->init_tile_id();
            app_node->set_mapped_to( t_id); // node knows to what tile;
            _tiles[ t_id].set_mapping( app_node);
        }
    }


    // (2) entertain user; compute overall results with std_dev.;
    long single_12_failures_count = _sketch_avg_migration_distance.size();
    double avg_migration_distance = 0;
    double avg_migration_distance_STDDEV = 0;   
    double avg_comm_volume_change = 0;
    double avg_comm_volume_change_STDDEV = 0;   
    for ( long i = 0; i < single_12_failures_count; i++) {
        avg_migration_distance += _sketch_avg_migration_distance[ i];
        avg_comm_volume_change += _sketch_avg_comm_volume_change[ i];
        printf ("\n %d migr_dist: %.2f comm_vol_change: %.2f", i,
            _sketch_avg_migration_distance[ i], _sketch_avg_comm_volume_change[ i]);
    }
    avg_migration_distance /= single_12_failures_count;
    avg_comm_volume_change /= single_12_failures_count;
    for ( long i = 0; i < single_12_failures_count; i++) {
        avg_migration_distance_STDDEV = avg_migration_distance_STDDEV +
            pow( (_sketch_avg_migration_distance[ i] - avg_migration_distance), 2);
        avg_comm_volume_change_STDDEV = avg_comm_volume_change_STDDEV +
            pow( (_sketch_avg_comm_volume_change[ i] - avg_comm_volume_change), 2);
    }
    avg_migration_distance_STDDEV /= single_12_failures_count;
    avg_comm_volume_change_STDDEV /= single_12_failures_count;
    avg_migration_distance_STDDEV = sqrt(avg_migration_distance_STDDEV);
    avg_comm_volume_change_STDDEV = sqrt(avg_comm_volume_change_STDDEV);

    printf ("\nSINGLE FAILURE INJECTION RESULTS");
    printf ("\n Avg. total migration distance: %.2f (STDDEV: %.2f)",
            avg_migration_distance, avg_migration_distance_STDDEV);
    printf ("\n Avg. comm volume change: %.2f (STDDEV: %.2f)\n",
            avg_comm_volume_change, avg_comm_volume_change_STDDEV);
}

bool REMNOC::run_remapping_two_simultaneous_failures_multiple_times()
{
    // same as run_remapping_single_failure_multiple_times(), but
    // two failures are injected at the same time;
    char msg[BUFFER_SIZE];
    double comm_volume_init = compute_total_comm_volume();
    _sketch_avg_migration_distance.clear();
    _sketch_avg_comm_volume_change.clear();

    // (1)
    long nodes_count = _app_graph->nodes_count();
    for ( long i1 = 0; i1 < nodes_count; i1++) {
        // () GUI stuff: entertain user before;
        if ( use_gui()) {
            sprintf( msg, "Initial mapping");
            gui()->update_screen( PRIORITY_MAJOR, msg, TILES);
        }       

        // () inject first failure into this node's tile/PE;
        _failures_count = 2;
        long t_id_failure_1 = _app_graph->get_node( i1)->tile_id();
        _tiles[ t_id_failure_1].set_is_broken( true);

        // inject the second failure as a simultaneous failure;
        for ( long i2 = i1 + 1; i2 < nodes_count; i2++) {
            // inject second failure;
            long t_id_failure_2 = _app_graph->get_node( i2)->tile_id();
            _tiles[ t_id_failure_2].set_is_broken( true);

            // () get the new mapping shape/polygon region (ideally convex);
            //run_computation_of_new_map_polygon_1();
            run_computation_of_new_map_polygon_2();
            // () assign cores to tiles inside the new mapping "shape";
            run_assignment_of_cores_to_new_locations();

            // () GUI stuff: entertain user after;
            if ( use_gui()) {
                sprintf( msg, "Mapping after second failure injection %d %d",i1,i2);
                gui()->update_screen( PRIORITY_MAJOR, msg, TILES);
            }

            // record the result of this i1 i2 concurrent failures;
            double total_migration_distance = double(compute_total_migration());
            double comm_volume_final = compute_total_comm_volume();
            double total_comm_volume_change = 100 * double(comm_volume_final - comm_volume_init)/
                double(comm_volume_init);
            _sketch_avg_migration_distance.push_back( total_migration_distance);
            _sketch_avg_comm_volume_change.push_back( total_comm_volume_change);

            // restore initial mapping we had before i1, i2 failures injection;
            _sketch_mapped_tiles_new.clear();
            _sketch_mapped_tiles_new.resize( nodes_count);
            for ( long sketch_i = 0; sketch_i < nodes_count; sketch_i++) {
                _sketch_mapped_tiles_new[ sketch_i] = 0;
            }
            for ( long j = 0; j < _tiles_count; j++) {
                _tiles[ j].set_mapping( 0);
            }
            // note I am backing up only i2 failure;
            _tiles[ t_id_failure_2].set_is_broken( false);
            for ( long j = 0; j < nodes_count; j++) {
                APPLICATION_NODE *app_node = _app_graph->get_node( j);
                long t_id = app_node->init_tile_id();
                app_node->set_mapped_to( t_id); // node knows to what tile;
                _tiles[ t_id].set_mapping( app_node);
            }
        }

        // () restore initial mapping to prepare for next i1 first failure injection;
        // _rects_black, _rects_white, _black_area, _white_area are reset
        // by run_computation_of_new_map_polygon_2();
        _sketch_mapped_tiles_new.clear();
        _sketch_mapped_tiles_new.resize( nodes_count);
        for ( long sketch_i = 0; sketch_i < nodes_count; sketch_i++) {
            _sketch_mapped_tiles_new[ sketch_i] = 0;
        }
        for ( long j = 0; j < _tiles_count; j++) {
            _tiles[ j].set_mapping( 0);
            _tiles[ j].set_is_broken( false);
        }   
        for ( long j = 0; j < nodes_count; j++) {
            APPLICATION_NODE *app_node = _app_graph->get_node( j);
            long t_id = app_node->init_tile_id();
            app_node->set_mapped_to( t_id); // node knows to what tile;
            _tiles[ t_id].set_mapping( app_node);
        }
    }

    // (2) entertain user; compute overall results with std_dev.;
    long single_12_failures_count = _sketch_avg_migration_distance.size();
    double avg_migration_distance = 0;
    double avg_migration_distance_STDDEV = 0;   
    double avg_comm_volume_change = 0;
    double avg_comm_volume_change_STDDEV = 0;   
    for ( long i = 0; i < single_12_failures_count; i++) {
        avg_migration_distance += _sketch_avg_migration_distance[ i];
        avg_comm_volume_change += _sketch_avg_comm_volume_change[ i];
        printf ("\n %d migr_dist: %.2f comm_vol_change: %.2f", i,
            _sketch_avg_migration_distance[ i], _sketch_avg_comm_volume_change[ i]);
    }
    avg_migration_distance /= single_12_failures_count;
    avg_comm_volume_change /= single_12_failures_count;
    for ( long i = 0; i < single_12_failures_count; i++) {
        avg_migration_distance_STDDEV = avg_migration_distance_STDDEV +
            pow( (_sketch_avg_migration_distance[ i] - avg_migration_distance), 2);
        avg_comm_volume_change_STDDEV = avg_comm_volume_change_STDDEV +
            pow( (_sketch_avg_comm_volume_change[ i] - avg_comm_volume_change), 2);
    }
    avg_migration_distance_STDDEV /= single_12_failures_count;
    avg_comm_volume_change_STDDEV /= single_12_failures_count;
    avg_migration_distance_STDDEV = sqrt(avg_migration_distance_STDDEV);
    avg_comm_volume_change_STDDEV = sqrt(avg_comm_volume_change_STDDEV);

    printf ("\nSINGLE FAILURE INJECTION RESULTS");
    printf ("\n Avg. total migration distance: %.2f (STDDEV: %.2f)",
            avg_migration_distance, avg_migration_distance_STDDEV);
    printf ("\n Avg. comm volume change: %.2f (STDDEV: %.2f)\n",
            avg_comm_volume_change, avg_comm_volume_change_STDDEV);
}

bool REMNOC::run_remapping_three_simultaneous_failures_multiple_times()
{
    // same as run_remapping_single_failure_multiple_times(), but
    // three failures are injected at the same time;
    char msg[BUFFER_SIZE];
    double comm_volume_init = compute_total_comm_volume();
    _sketch_avg_migration_distance.clear();
    _sketch_avg_comm_volume_change.clear();

    // (1) inject three failures multiple times;
    long nodes_count = _app_graph->nodes_count();
    for ( long i1 = 0; i1 < nodes_count; i1++) {
        // () GUI stuff: entertain user before;
        if ( use_gui()) {
            sprintf( msg, "Initial mapping");
            gui()->update_screen( PRIORITY_MAJOR, msg, TILES);
        }       

        // () inject first failure into this node's tile/PE;
        _failures_count = 3;
        long t_id_failure_1 = _app_graph->get_node( i1)->tile_id();
        _tiles[ t_id_failure_1].set_is_broken( true);

        // inject second failure as a simultaneous failure;
        for ( long i2 = i1 + 1; i2 < nodes_count; i2++) {
            // inject 2nd failure;
            long t_id_failure_2 = _app_graph->get_node( i2)->tile_id();
            _tiles[ t_id_failure_2].set_is_broken( true);

            // inject third failure as a simultaneous failure;
            for ( long i3 = i2 + 1; i3 < nodes_count; i3++) {
                // inject 3rd failure;
                long t_id_failure_3 = _app_graph->get_node( i3)->tile_id();
                _tiles[ t_id_failure_3].set_is_broken( true);

                // () get the new mapping shape/polygon region (ideally convex);
                //run_computation_of_new_map_polygon_1();
                run_computation_of_new_map_polygon_2();
                // () assign cores to tiles inside the new mapping "shape";
                run_assignment_of_cores_to_new_locations();

                // () GUI stuff: entertain user after;
                if ( use_gui()) {
                    sprintf( msg, "Mapping after three failures injection");
                    gui()->update_screen( PRIORITY_MAJOR, msg, TILES);
                }

                // record the result of this three concurrent failures;
                double total_migration_distance = double(compute_total_migration());
                double comm_volume_final = compute_total_comm_volume();
                double total_comm_volume_change = 100 * double(comm_volume_final - comm_volume_init)/
                    double(comm_volume_init);
                _sketch_avg_migration_distance.push_back( total_migration_distance);
                _sketch_avg_comm_volume_change.push_back( total_comm_volume_change);

                // restore initial mapping;
                _sketch_mapped_tiles_new.clear();
                _sketch_mapped_tiles_new.resize( nodes_count);
                for ( long sketch_i = 0; sketch_i < nodes_count; sketch_i++) {
                    _sketch_mapped_tiles_new[ sketch_i] = 0;
                }
                for ( long j = 0; j < _tiles_count; j++) {
                    _tiles[ j].set_mapping( 0); 
                }
                _tiles[ t_id_failure_3].set_is_broken( false);
                for ( long j = 0; j < nodes_count; j++) {
                    APPLICATION_NODE *app_node = _app_graph->get_node( j);
                    long t_id = app_node->init_tile_id();
                    app_node->set_mapped_to( t_id); // node knows to what tile;
                    _tiles[ t_id].set_mapping( app_node);
                }
            }
            _tiles[ t_id_failure_2].set_is_broken( false);
        }
        _tiles[ t_id_failure_1].set_is_broken( false);
    }

    // (2) entertain user; compute overall results with std_dev.;
    long multiple_3_failures_count = _sketch_avg_migration_distance.size();
    double avg_migration_distance = 0;
    double avg_migration_distance_STDDEV = 0;   
    double avg_comm_volume_change = 0;
    double avg_comm_volume_change_STDDEV = 0;   
    for ( long i = 0; i <  multiple_3_failures_count; i++) {
        avg_migration_distance += _sketch_avg_migration_distance[ i];
        avg_comm_volume_change += _sketch_avg_comm_volume_change[ i];
        printf ("\n %d migr_dist: %.2f comm_vol_change: %.2f", i,
            _sketch_avg_migration_distance[ i], _sketch_avg_comm_volume_change[ i]);
    }
    avg_migration_distance /=  multiple_3_failures_count;
    avg_comm_volume_change /=  multiple_3_failures_count;
    for ( long i = 0; i <  multiple_3_failures_count; i++) {
        avg_migration_distance_STDDEV = avg_migration_distance_STDDEV +
            pow( (_sketch_avg_migration_distance[ i] - avg_migration_distance), 2);
        avg_comm_volume_change_STDDEV = avg_comm_volume_change_STDDEV +
            pow( (_sketch_avg_comm_volume_change[ i] - avg_comm_volume_change), 2);
    }
    avg_migration_distance_STDDEV /=  multiple_3_failures_count;
    avg_comm_volume_change_STDDEV /=  multiple_3_failures_count;
    avg_migration_distance_STDDEV = sqrt(avg_migration_distance_STDDEV);
    avg_comm_volume_change_STDDEV = sqrt(avg_comm_volume_change_STDDEV);

    printf ("\nSINGLE FAILURE INJECTION RESULTS");
    printf ("\n Avg. total migration distance: %.2f (STDDEV: %.2f)",
            avg_migration_distance, avg_migration_distance_STDDEV);
    printf ("\n Avg. comm volume change: %.2f (STDDEV: %.2f)\n",
            avg_comm_volume_change, avg_comm_volume_change_STDDEV);
}

void REMNOC::run_computation_of_new_map_polygon_1()
{
    // get the new mapping polygon and put in _sketch_mapped_tiles_new;
    // idea: do simultaneous ver and hor stripifying using coordinates of 
    // chip, initial mapping and broken tiles;
    // black rects are those that are good and white ones are the empty
    // ones; they can be of different areas; re-mapping to find map-shape 
    // means converting some whites to black and prunning away some blacks;


    // (a) build the rects used for phase one: finding new region with
    // good tiles only, as close as possible to the initial one and with
    // max overlap for minimum disruption;
    build_rects_arrays();

    
    // (b)
    compute_distances_white_to_black();
    sort( _rects_black.begin(), _rects_black.end(), compare_by_area);
    sort( _rects_white.begin(), _rects_white.end(), compare_by_dist);
    // by this time _rects_white contains rects that are smalest
    // distance to any black and in decreasing order of their size;
    // from end to front; example of (size,dist) _rects_white:
    // {(1,3)(1,3)(2,3)  (5,2)  (1,1)(3,1)(3,1)}
    // _rects_white will be moved to blacks:(3,1),(3,1)...
    // _rects_black are in nonincreasing order of their area;
    print_sketch_arrays(); // debug;


    // (c)
    _app_area = _app_graph->nodes_count();
    

    // (d) move rect from white set to black set to make up for
    // displaced area;
    _black_area = _app_area - _displaced_area;
    printf("\n _app_area=%d _black_area=%d", _app_area, _black_area);
    vector<RECT> buffer_rects_black;
    while ( _black_area < _app_area) {
        RECT w_rect = _rects_white[ _rects_white.size() - 1];
        printf(" w_rect_area=%d", w_rect.area());
        _black_area += w_rect.area();
        _white_area -= w_rect.area();
        buffer_rects_black.push_back( w_rect); // put in blacks buffer;
        _rects_white.pop_back(); // remove from whites;
    }
    printf("\n _app_area=%d _black_area=%d", _app_area, _black_area);
    // if now the black area is more than is needed, trim some of the
    // appendices with area only 1, for now; later we may want to remove
    // appendices with larger area;
    int appendix_area_threshold = 1;
    if ( _black_area > _app_area) {
        int surplus_area = _black_area - _app_area;
        while ( surplus_area > 0) {
            if ( _rects_black[ _rects_black.size() - 1].area() <= appendix_area_threshold) {
                RECT b_rect = _rects_black[ _rects_black.size() - 1];
                _black_area -= b_rect.area();
                _white_area += b_rect.area();
                _rects_white.push_back( b_rect);
                _rects_black.pop_back(); // remove from blacks;
            }
            int old_surplus_area = surplus_area;
            surplus_area = _black_area - _app_area;
            if (surplus_area == old_surplus_area) break;
        }
    }
    for ( int i = 0; i < buffer_rects_black.size(); i++) {
        _rects_black.push_back( buffer_rects_black[ i]);
    }


    // (e) by now the black rects that make-up the new mapping-shape
    // are computed; have to store the "covered" good tiles id's into
    // _sketch_mapped_tiles_new array;
    int nodes_count = int(_app_graph->nodes_count());
    int black_rects_count = _rects_black.size();
    int counter = 0;
    for ( int i = 0; i < black_rects_count; i ++) {
        int xl = _rects_black[ i].xl();
        int xr = _rects_black[ i].xr();     
        int yb = _rects_black[ i].yb();
        int yt = _rects_black[ i].yt();
        for ( int y = yb; y < yt; y++) {
            for ( int x = xl; x < xr; x++) {
                int t_id = x + y * _nx;
                if ( counter < nodes_count) {
                    _sketch_mapped_tiles_new[ counter] = t_id;
                } else {
                    // increase a little the new mapping shape if more
                    // tiles "came in" due to white rects moved to blacks;
                    _sketch_mapped_tiles_new.push_back(t_id);                   
                }
                counter ++;
            }
        }
    }
}

void REMNOC::run_computation_of_new_map_polygon_2()
{
    // get the new mapping polygon and put in _sketch_mapped_tiles_new;
    // idea: use center of gravity of all current black rects; convert 
    // to blacks the closest white to the center of gravity of the blacks;
    // do it for |nodes_count / 2|;
    // Warning: this has longer runtime?


    // ()
    build_rects_arrays_with_unit_area(); // blacks and whites;


    // ()
    _app_area = _app_graph->nodes_count();


    // ()
    int converted_count = 0;
    while ( converted_count < _failures_count) {
        converted_count ++;
        // update center of mass of blacks;
        //compute_mass_center_black(); // simply mass center;
        compute_mass_center_black_weighted(); // weighted rects mass center;
        
        // update distance of white to new center of mass;
        compute_Euclidean_distances_white_to_mass_center_black();
        // get last one with smallest distance from the blacks center of mass;
        sort( _rects_white.begin(), _rects_white.end(), compare_by_dist);
        //print_sketch_arrays(); // debug;
        // convert the closest white to black;
        RECT w_rect = _rects_white[ _rects_white.size() - 1];
        _black_area += w_rect.area();
        _white_area -= w_rect.area();
        _rects_black.push_back( w_rect); // put in blacks buffer;
        _rects_white.pop_back(); // remove from whites;
    }

    // () by now the black rects that make-up the new mapping-shape
    // are computed; have to store the "covered" good tiles id's into
    // _sketch_mapped_tiles_new array;
    int nodes_count = int(_app_graph->nodes_count());
    int black_rects_count = _rects_black.size();
    int counter = 0;
    for ( int i = 0; i < black_rects_count; i ++) {
        int xl = _rects_black[ i].xl();
        int xr = _rects_black[ i].xr();     
        int yb = _rects_black[ i].yb();
        int yt = _rects_black[ i].yt();
        for ( int y = yb; y < yt; y++) {
            for ( int x = xl; x < xr; x++) {
                int t_id = x + y * _nx;
                if ( counter < nodes_count) {
                    _sketch_mapped_tiles_new[ counter] = t_id;
                } else {
                    // increase a little the new mapping shape if more
                    // tiles "came in" due to white rects moved to blacks;
                    _sketch_mapped_tiles_new.push_back(t_id);                   
                }
                counter ++;
            }
        }
    }
}

void REMNOC::compute_distances_white_to_black()
{
    int b_count = _rects_black.size();
    int w_count = _rects_white.size();
    for ( int j = 0; j < w_count; j++) {
        int min_dist = INT_MAX; // of this white guy;
        for ( int i = 0; i < b_count; i++) {
            int this_dist = my_white_black_distance( j, i);
            if ( this_dist < min_dist) {
                min_dist = this_dist;
            }
            if ( min_dist == 0) break;
        }
        _rects_white[ j].set_dist( double(min_dist));
    }
}

void REMNOC::compute_Euclidean_distances_white_to_mass_center_black()
{
    // compute E distance of each white rect to the mass center of 
    // the black rects;
    int w_count = _rects_white.size();
    for ( int j = 0; j < w_count; j++) {
        double this_dist = my_Euclidean_white_mass_center_black_distance( j);
        _rects_white[ j].set_dist( this_dist);
    }
}

double REMNOC::my_Euclidean_white_mass_center_black_distance( int j)
{
    // j is index of rect in _rects_white;
    double x_center_w = double(_rects_white[ j].xl() + _rects_white[ j].xr()) / 2;
    double y_center_w = double(_rects_white[ j].yb() + _rects_white[ j].yt()) / 2;
    double result = 
        sqrt( (x_center_w - _x_mass_center_black)*(x_center_w - _x_mass_center_black) +
            (y_center_w - _y_mass_center_black)*(y_center_w - _y_mass_center_black));
    return (result);
}

int REMNOC::my_white_black_distance( int j, int i)
{
    // j is index of rect in _rects_white and i is index of 
    // rect in _rects_black; r2 (white) can be in any of the 
    // following 8 regions, relativ to r1 (black):
    // Region 3 | R4 |  R5
    // ---------|----|-----
    // R2       | r1 |  R6
    // ---------|----|-----
    // R1       | R8 |  R7
    int result = 0;
    int xl_1 = _rects_black[ i].xl();
    int xr_1 = _rects_black[ i].xr();
    int yb_1 = _rects_black[ i].yb();
    int yt_1 = _rects_black[ i].yt();
    int xl_2 = _rects_white[ j].xl();
    int xr_2 = _rects_white[ j].xr();
    int yb_2 = _rects_white[ j].yb();
    int yt_2 = _rects_white[ j].yt();

    float delta_x = 0, delta_y = 0, delta = 0;
    if ( xl_2 < xl_1) {
        if ( yb_2 < yb_1) { // R1
            delta_x = (xr_2 < xl_1) ? (xl_1 - xr_2) : 0;
            delta_y = (yt_2 < yb_1) ? (yb_1 - yt_2) : 0;
            delta = delta_x + delta_y + 1;
            // +1 to encourage relocation to West or South rather than SW;
        } else if ( yb_2 >= yb_1 && yb_2 < yt_1) { // R2
            delta_x = xl_1 - xr_2;
            delta_y = 0;
            delta = delta_x + delta_y;
        } else if ( yb_2 >= yt_1) { // R3
            delta_x = (xr_2 < xl_1) ? (xl_1 - xr_2) : 0;
            delta_y = yb_2 - yt_1;
            delta = delta_x + delta_y + 1;
        }
    } else if ( xl_2 >= xl_1 && xl_2 < xr_1) {
        if ( yt_2 <= yb_1) { // R8
            delta_x = 0;
            delta_y = yb_1 - yt_2;
            delta = delta_x + delta_y;
        } else if ( yb_2 >= yt_1) { // R4
            delta_x = 0;
            delta_y = yb_2 - yt_1;
            delta = delta_x + delta_y;
        }
    } else if ( xl_2 >= xr_1) {
        if ( yb_2 < yb_1) { // R7
            delta_x = (xl_2 > xr_1) ? (xl_2 - xr_1) : 0;
            delta_y = (yt_2 < yb_1) ? (yb_1 - yt_2) : 0;
            delta = delta_x + delta_y + 1;
        } else if ( yb_2 >= yb_1 && yb_2 < yt_1) { // R6
            delta_x = xl_2 - xr_1;
            delta_y = 0;
            delta = delta_x + delta_y;
        } else if ( yb_2 >= yt_1) { // R5
            delta_x = (xl_2 > xr_1) ? (xl_2 - xr_1) : 0;
            delta_y = yb_2 - yt_1;
            delta = delta_x + delta_y + 1;
        }
    }
    result = int(delta);
    return result;
}

void REMNOC::run_assignment_of_cores_to_new_locations()
{
    // _sketch_mapped_tiles_new now should contain tiles (not all tiles) id's which 
    // form the new mapping-shape, region where cores have to be mapped;
    // this region should be as convex as possible and excluding broken tiles;

    long nodes_count = _app_graph->nodes_count();
    long mapped_tiles_new_count = _sketch_mapped_tiles_new.size(); // may be > nodes_count;
    long old_x = 0, old_y = 0;
    long new_x = 0, new_y = 0;
    long cost = 0;

    HUNGARIAN_ONE hungarian; // will solve the linear assignment problem;
    // cores count, mapped-to-tiles count = cores count, but we could consider
    // all tiles of the NOC, and penalize more assignment of cores to tiles
    // "outside" of the new mapping-shape/boundaries;
    hungarian.initialize( nodes_count, mapped_tiles_new_count); 

    for ( long i = 0; i < nodes_count; i++) {
        APPLICATION_NODE *app_node = _app_graph->get_node( i);
        long old_tile_id = app_node->tile_id();
        old_x = _tiles[ old_tile_id].x();
        old_y = _tiles[ old_tile_id].y();

        // compute the cost of assigning this core to any of the 
        // new locations, within the shape of the new mapping;
        for ( long j = 0; j < mapped_tiles_new_count; j++) {
            new_x = _tiles[ _sketch_mapped_tiles_new[ j]].x();
            new_y = _tiles[ _sketch_mapped_tiles_new[ j]].y();
            long dist = abs(new_x - old_x) + abs(new_y - old_y);
            //cost = dist * dist * TASK_MIGRATION_ENERGY_COST_PER_HOP;
            cost = dist * 1 * long(app_node->io_comm_volume());

            hungarian.set_cost(i, j, cost);
        }
    }
    // call the actual magyar man;
    hungarian.run_hungarian();

    // clean previous mapping in nodes and tiles;
    clean_tiles_mappings();

    // get new assignments and update nodes/cores and tiles they are mapped to;
    for ( long i = 0; i < nodes_count; i ++) {
        long new_tile_id = 
            _sketch_mapped_tiles_new[ hungarian.report_assignment_of( i)];

        APPLICATION_NODE *app_node = _app_graph->get_node( i);
        _tiles[ new_tile_id].set_mapping( app_node); // update tile's ptr to core;
        app_node->set_mapped_to( new_tile_id);
    }   

    //hungarian.print_hungarian_assignment(); // debug;
}

long REMNOC::compute_total_migration()
{
    // of all cores, as aggregated total Manhatan distance from old to new
    // locations; return total number of hops IP/cores will have to be
    // shuffled;
    long result = 0;
    long old_x = 0, old_y = 0;
    long new_x = 0, new_y = 0;
    long nodes_count = _app_graph->nodes_count();
    for ( long i = 0; i < nodes_count; i++) {
        APPLICATION_NODE *node = _app_graph->get_node( i);
        long old_tile_id = node->init_tile_id();
        old_x = _tiles[ old_tile_id].x();
        old_y = _tiles[ old_tile_id].y();
        long new_tile_id = node->tile_id();
        new_x = _tiles[ new_tile_id].x();
        new_y = _tiles[ new_tile_id].y();
        result += abs(new_x - old_x) + abs(new_y - old_y);
    }
    return result;
}

double REMNOC::compute_total_comm_volume()
{
    // of all arcs; return aggregated comm volume weighted by manhattan 
    // length of each arc/link;
    double result = 0;
    long src_x = 0, src_y = 0, des_x = 0, des_y = 0;
    long arcs_count = _app_graph->arcs_count();
    for ( long j = 0; j < arcs_count; j++) {
        APPLICATION_ARC *arc = _app_graph->get_arc( j);
        APPLICATION_NODE *src_node = _app_graph->get_node( arc->src_id());
        APPLICATION_NODE *des_node = _app_graph->get_node( arc->des_id());
        long src_tile_id = src_node->tile_id();
        src_x = _tiles[ src_tile_id].x();
        src_y = _tiles[ src_tile_id].y();
        long des_tile_id = des_node->tile_id();
        des_x = _tiles[ des_tile_id].x();
        des_y = _tiles[ des_tile_id].y();
        long md = abs(des_x - src_x) + abs(des_y - src_y);
        result += md * arc->comm_volume();
    }
    return result;
}

double REMNOC::compute_total_comm_volume_for_annealer( int move_type)
{
    // of all arcs; return aggregated comm volume weighted by manhattan 
    // length of each arc/link "+" distance from preferred location;
    // move_type: "0" move, "1" swap; for swap mode, then I do not penalize 
    // new location far away from initial location in the initial mapping,
    // because in the second annealing, refining process, I move nodes
    // within the new mapping region in order to improve total comm volume,
    // displacement is not penalized;

    double result1 = 0;
    long src_x = 0, src_y = 0, des_x = 0, des_y = 0;
    long arcs_count = _app_graph->arcs_count();
    for ( long j = 0; j < arcs_count; j++) {
        APPLICATION_ARC *arc = _app_graph->get_arc( j);
        APPLICATION_NODE *src_node = _app_graph->get_node( arc->src_id());
        APPLICATION_NODE *des_node = _app_graph->get_node( arc->des_id());
        long src_tile_id = src_node->tile_id();
        src_x = _tiles[ src_tile_id].x();
        src_y = _tiles[ src_tile_id].y();
        long des_tile_id = des_node->tile_id();
        des_x = _tiles[ des_tile_id].x();
        des_y = _tiles[ des_tile_id].y();
        long md = abs(des_x - src_x) + abs(des_y - src_y);
        result1 += md * arc->comm_volume();
    }

    // add also the penalty due to "distance" from initial mapping;
    double result2 = 0;
    if ( move_type == 0) { // migration/displacement is penalized;
        long old_x = 0, old_y = 0, new_x = 0, new_y = 0;
        long nodes_count = _app_graph->nodes_count();
        for ( long i = 0; i < nodes_count; i++) {
            APPLICATION_NODE *node = _app_graph->get_node( i);
            long old_tile_id = node->init_tile_id();
            old_x = _tiles[ old_tile_id].x();
            old_y = _tiles[ old_tile_id].y();
            long new_tile_id = node->tile_id();
            new_x = _tiles[ new_tile_id].x();
            new_y = _tiles[ new_tile_id].y();
            double dist = sqrt( (new_x - old_x)*(new_x - old_x) +
                                (new_y - old_y)*(new_y - old_y)); // Euclidean dist;
            result2 += dist * dist *
                node->io_comm_volume(); // weight by outdegree comm volume of node;
        }
    }
    
    return (result1 + result2);
}

void REMNOC::check_remapping()
{
}

////////////////////////////////////////////////////////////////////////////////
//
// APPLICATION_GRAPH
//
////////////////////////////////////////////////////////////////////////////////

void APPLICATION_GRAPH::print_graph()
{
    // debug;
    printf("\nApplication graph %d", _id);
    printf("\nIP/cores: %d", _nodes.size());
    for ( long i = 0; i < _nodes_count; i ++) {
        printf("\n%d\t%s", i, _nodes[i].name().c_str());
        printf("        \tfin:");
        for ( long k = 0; k < _nodes[i].fanin().size(); k ++) {
            printf(" %d", _nodes[i].fanin()[k]);
        }
        printf("\tfout:");
        for ( long k = 0; k < _nodes[i].fanout().size(); k ++) {
            printf(" %d", _nodes[i].fanout()[k]);
        }
        printf("\tcomm_vol: %.1f", _nodes[i].io_comm_volume());
    }
    printf("\nArcs: %d", _arcs.size());
    for ( long j = 0; j < _arcs_count; j ++) {
        printf("\n%d -> %d    \t%d", 
               _arcs[j].src_id(), _arcs[j].des_id(), _arcs[j].comm_volume());
    }
    printf("\n");
}
