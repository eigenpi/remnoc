////////////////////////////////////////////////////////////////////////////////
//
// Author:     Cristinel Ababei
// E-mail:     cristinel.ababei@ndsu.edu
// Copyright:  Please see README.txt
//
////////////////////////////////////////////////////////////////////////////////

#include "remnoc_hungarian.h"
#include "remnoc_gui.h"
#include "remnoc.h"
#include <assert.h>
#include <time.h>

using namespace std;

////////////////////////////////////////////////////////////////////////////////
//
// launching point;
//
////////////////////////////////////////////////////////////////////////////////

int main( int argc, char *argv[])
{
    // welcome;
    char welcome[] =
        "\n--------------------------------------------------\n"
        "REMNOC: REMapping NOC\n"
        "cristinel.ababei@ndsu.edu, Compiled "__DATE__" \n"
        "--------------------------------------------------\n";
    printf("%s", welcome);

    time_t start_sec, end_sec; // seconds;
    timeval time_val_1, time_val_2; // sec and microsec;
    clock_t start_clock, end_clock;
    clock_t diff_clock;
    time_t diff_sec;    
    char msg[BUFFER_SIZE];
 

    // (1) create the network and then readin/record the current mapping
    // from mapping_file.map;
    REMNOC remnoc;

    
    // (2) create the application graph(s) from the app_file.app; for now
    // we work with one application only;
    remnoc.parse_command_arguments( argc, argv);
    remnoc.create_application_graphs( argc, argv);
    // read current mapping(s) of app graph(s);
    remnoc.read_in_mappings( argc, argv);
    remnoc.print_initial_stats( argc, argv); // entertain user;
    remnoc.print_network_tiles(); // debug;
    double comm_volume_init = remnoc.compute_total_comm_volume();
    
    // (3) create empty gui object and populate if required by user;
    GUI_GRAPHICS gui( &remnoc);
    remnoc.set_gui( &gui); // innitially gui is empty;
    if ( remnoc.use_gui()) { // if user asked to use the gui;
        // mark flag that we are gonna use the gui; set is_gui_usable
        // and wait_for_user_input_automode; then build;
        gui.set_graphics_state( true, 1);
        gui.build();
        gui.init_draw_coords( 100.0, 100.0);
    } else { // gui is not usable;
        gui.set_graphics_state( false, 1);
    }



    bool regular_run = true;
    // regular run is when command arguments are used and the remapping
    // is run once only; for collecting results, we inject single or
    // multiple failures for all possible combinations and then get 
    // avg and stddev.; collecting results happens in the "else" branch;
    if ( regular_run) {

        // (4) GUI stuff: entertain user before;
        if ( remnoc.use_gui()) {
            sprintf( msg, "Initial mappings %d", 1);
            remnoc.gui()->update_screen( PRIORITY_MAJOR, msg, TILES);
        }

        // (5) injects failures into selected tiles;
        remnoc.inject_failures();

        // (6) do remapping and measure time spend only on it;
        start_sec = time(NULL);
        gettimeofday( &time_val_1, 0);
        start_clock = clock();
        assert(start_sec != (time_t)(-1));
        assert(start_clock != (clock_t)(-1));

        // run the thing;
        if ( remnoc.remapping_algo() == REMAPPING_HEURISTIC) {
            remnoc.run_remapping();
        } else {
            remnoc.run_simulated_annealing_remapping();
        }
        
        end_sec  = time(NULL);
        gettimeofday( &time_val_2, 0);
        end_clock = clock();
        assert(end_sec != (time_t)(-1));
        assert(end_clock != (clock_t)(-1));
        diff_sec = end_sec - start_sec;
        diff_clock = end_clock - start_clock;


        // (7) GUI stuff: entertain user after;
        if ( remnoc.use_gui()) {
            sprintf( msg, "Final mappings %d", 1);
            remnoc.gui()->update_screen( PRIORITY_MAJOR, msg, TILES);
        }

        // (8) entertain user - text mode;
        remnoc.print_network_tiles(); // debug;
        long total_migration_distance = remnoc.compute_total_migration();   
        double comm_volume_final = remnoc.compute_total_comm_volume();
        float total_comm_volume_change = 100 * float(comm_volume_final - comm_volume_init)/
            float(comm_volume_init);

        printf ("\n");
        //printf ("cputime : start_clock = %lu units, end_clock = %lu units\n", start_clock, end_clock);
        printf ("cputime : processor time used = %2.4f sec\n", (double)diff_clock/CLOCKS_PER_SEC);
        //printf ("walltime : start_time = %lu sec, end_time = %lu sec\n", start_sec, end_sec);
        //printf ("walltime : elapsed (wall clock) time = %lu sec\n", diff_sec);
        double diff_sec_usec = time_val_2.tv_sec - time_val_1.tv_sec + 
            double(time_val_2.tv_usec - time_val_1.tv_usec) / 1000000.0;
        printf ("walltime : elapsed (wall clock) time = %2.8f sec\n", diff_sec_usec);
        printf ("\n");
        printf ("Total migration distance : %d\n", total_migration_distance);
        printf ("Comm volume : init %.1f  final %.1f  change %.2f \%\n",
                comm_volume_init, comm_volume_final, total_comm_volume_change);
    }



    // "-mode" is to collect results;
    else {
        // (1) run the thing multiple times using proposed algo;
        // gui stuff is done inside this call;
        //remnoc.run_remapping_single_failure_multiple_times();
        //remnoc.run_remapping_two_sequential_failures_multiple_times();
        //remnoc.run_remapping_two_simultaneous_failures_multiple_times();
        //remnoc.run_remapping_three_simultaneous_failures_multiple_times();

        // (2) run the thing multiple times using the simulated annealing;
        remnoc.run_remapping_single_failure_multiple_times_SA_based();

    }


        
    // (9) do some clean-up;
    if ( gui.is_gui_usable()) {
        gui.close_graphics(); // close down X Display;
    }
}

